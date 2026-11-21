// Reference for the embedded team (ASTRA abort integration, TOAD/cpp/INTEGRATION.md); not compiled in this repo.
// TrajectoryFollower.cpp from the FC firmware, modified to fly the PolyQP abort at the "autoland" branch. Only the lines
// marked "ABORT:" changed; everything else is the original follower.
#include "TrajectoryFollower.h"
#include "ASTRA_Trajectory_Data.h"
#include "Abort.h" // ABORT: the abort harness (the only abort header the FC includes)
#include "CommandRouter.h"
#include "CommsSerial.h"
#include "FlashLogging.h"
#include "GPS.h"
#include "GimbalServos.h"
#include "IMU.h"
#include "Mag.h"
#include "Prop.h"
#include "RollControl.h"
#include "TrajectoryLogger.h"
#include "astra_structs.h"
#include "controller_and_estimator.h"
#include "elapsedMillis.h"
#include "fc_pins.h"
#include <Arduino.h>

#define TELEMETRY_INTERVAL_US 100000
#define COMMAND_INTERVAL_US 1000
#define LOG_X_EST_INTERVAL_US 100000
#define HB_KILL_INTERVAL_MS 500 // if we don't see a heartbeat command every 0.5 seconds during flight, disable gimbal and props

namespace TrajectoryFollower {

// ABORT: option (b) of INTEGRATION.md: the nominal Trajectory::TRAJ_* arrays stay const in flash, the abort is written into
// this one static RAM buffer (MAX_POINTS x 80 bytes = 48 KB), and the follower reads through pointers that switch on abort.
static float ABORT_TIME[Abort::MAX_POINTS];
static float ABORT_STATES[Abort::MAX_POINTS][15];
static float ABORT_INPUTS[Abort::MAX_POINTS][4];

bool kill_flag;
bool arm_flag;
bool auto_land_cmd;
unsigned long last_hb;
elapsedMicros timer;

float lerp(float a, float b, float t0, float t1, float t) {
  if (t <= t0)
    return a;
  if (t >= t1)
    return b;
  if (t0 == t1)
    return b; // immediately get to b
  return a + (b - a) * ((t - t0) / (t1 - t0));
}

void follow_trajectory() {
  bool flight_armed = false;
  kill_flag = false;
  arm_flag = false;
  auto_land_cmd = false;

  long counter = 0;

  GimbalServos::centerGimbal();
  Mag::beginMeasurement();

  GPS_Point last_gps_pos = {-1, -1, -1}; // first packet will be marked as new
  ControllerAndEstimator::init_controller_and_estimator_constants();
  Controller_Internals cs;
  telemetry_packet_t fp;

  while (!Mag::isMeasurementReady()) {
    CommsSerial.println("Waiting on mag...");
    delay(100);
  }
  while (!GPS::has_valid_recent_pos()) {
    GPS::pump_events();
    CommsSerial.println("Waiting on gps...");
    delay(100);
  }

  GPS::set_current_position_as_origin();

  TrajectoryLogger::log_calib_flash();

  timer = elapsedMicros();
  unsigned long lasttelemetry = timer;
  unsigned long lastloop = timer;
  unsigned long lastlogx_est = timer;
  unsigned long led_on_time = 0;
  bool led_on = false;
  bool auto_land = false;
  last_hb = timer;

  float last_time_s = timer / 1000000.0;
  float mx, my, mz;

  // ABORT: read pointers into the trajectory being flown (the nominal one until an abort switches them to the abort buffer)
  const float *traj_time = Trajectory::TRAJ_TIME;
  const float (*traj_states)[15] = Trajectory::TRAJ_STATES;
  const float (*traj_inputs)[4] = Trajectory::TRAJ_INPUTS;
  int traj_num_points = Trajectory::TRAJ_NUM_POINTS;
  float traj_t0 = 0.0f;           // ABORT: time base of the trajectory being flown (s since arm; the abort re-bases it)
  Controller_Output last_co = {}; // ABORT: previous tick's command, for the abort's start-acceleration boundary condition
  int pt = 0;

  // ABORT: one loop that advances pt by time replaces the nested for (pt) / while, so pt and the time base can be reset
  // mid-flight. Same behaviour as before: pt stays 0 until armed, then moves to the segment containing the time of the
  // previous tick, and the loop ends once that time passes the last point of the trajectory being flown.
  while (true) {
    if (flight_armed) {
      while (pt < traj_num_points - 1 && last_time_s - traj_t0 >= traj_time[pt + 1]) {
        pt++;
      }
      if (pt >= traj_num_points - 1) {
        break;
      }
    }
    { // ABORT: the former inner while body follows, unchanged except where marked
      while (CommsSerial.available()) {
        CommandRouter::receive_byte(CommsSerial.read());
      }

      if (!auto_land & (auto_land_cmd || ((timer - last_hb) / 1000.0 > HB_KILL_INTERVAL_MS))) {
        if (flight_armed) {
          auto_land_cmd = false;
          CommsSerial.println("Autoland sequence activated!");
          // ABORT: plan from the estimate the controller used on the previous tick, the previous tick's command and the LESO
          // disturbance estimate. u_cmd follows LESO_Position / ASTRAv2_Controller: [U(0); U(1); thrust N; roll]. This
          // controller_and_estimator.cpp writes U(0) into co.gimbal_yaw_deg and U(1) into co.gimbal_pitch_deg (after its
          // gimbal EMA filter), so the two are taken in that order; check this against the FC copy (INTEGRATION.md).
          unsigned long plan_start_us = timer;
          Vector4 u_prev(last_co.gimbal_yaw_deg * DEG_TO_RAD, last_co.gimbal_pitch_deg * DEG_TO_RAD, last_co.thrust_N, last_co.roll_rad_sec_squared);
          Abort::Result abort_result = Abort::plan(ControllerAndEstimator::get_controller_x_est(), u_prev, ControllerAndEstimator::get_position_disturbance(),
                                                   ABORT_TIME, ABORT_STATES, ABORT_INPUTS, Abort::MAX_POINTS);
          // ABORT: fly the abort buffer from its t = 0, which is the previous tick (the estimate it planned from). The solve
          // blocks this loop, so the follower resumes partway into the plan, where the vehicle is by then.
          traj_time = ABORT_TIME;
          traj_states = ABORT_STATES;
          traj_inputs = ABORT_INPUTS;
          traj_num_points = abort_result.num_points;
          traj_t0 = last_time_s;
          pt = 0;
          lastloop = timer; // ABORT: the blocking solve overran this tick; restart the 1 kHz schedule instead of catching up
          CommsSerial.printf("Abort plan: status %u (0 OK, 1 INFEASIBLE_FLOWN, 2 FALLBACK), solver %u, %d points, %.2f s, solved in %.1f ms\n",
                             static_cast<unsigned>(abort_result.status), static_cast<unsigned>(abort_result.solver_status), abort_result.num_points, abort_result.duration_s,
                             (timer - plan_start_us) / 1000.0); // ABORT: log the Result
          auto_land = true;
        } else {
          CommsSerial.println("Radio disconnected pre-flight");
          kill_flag = true;
        }
      }

      if (kill_flag) {
        break;
      }

      if (arm_flag) {
        arm_flag = false; // only do this on rising edge
        flight_armed = true;
        timer = elapsedMicros();
        last_time_s = timer / 1000000.0;
        lasttelemetry = timer;
        lastloop = timer;
        lastlogx_est = timer;
        last_hb = timer;
        counter = 0;
        GPS::set_current_position_as_origin();

        TrajectoryLogger::log_x_est();
        TrajectoryLogger::log_flight_p(); // only log initial controller state, this is too much data to log every frame

        led_on_time = millis();
        led_on = true;
        digitalWrite(PIN_TEST_LED, LOW);
      }

      if (led_on && millis() - led_on_time > 500) {
        digitalWrite(PIN_TEST_LED, HIGH);
      }

      Controller_Input ci;

      IMU::Data imu_reading;
      IMU::IMUs[0].read_latest(&imu_reading);

      if (Mag::isMeasurementReady()) {
        ci.new_mag_packet = true;
        Mag::read_xyz_calibrated(mx, my, mz);
        Mag::beginMeasurement();
      } else {
        ci.new_mag_packet = false;
      }

      GPS::pump_events();
      GPS_Point gps_rel_pos = GPS::get_rel_xyz_pos();
      GPS_Velocity gps_vel = GPS::get_velocity();
      GPS::fill_covariances(&ci.gps);

      ci.imu.accel_x = -imu_reading.acc[2];
      ci.imu.accel_y = imu_reading.acc[1];
      ci.imu.accel_z = imu_reading.acc[0];
      ci.imu.gyro_yaw = -imu_reading.gyro[2];
      ci.imu.gyro_pitch = imu_reading.gyro[1];
      ci.imu.gyro_roll = imu_reading.gyro[0];
      ci.mag.mag_x = -mz;
      ci.mag.mag_y = my;
      ci.mag.mag_z = mx;
      ci.gps.pos = gps_rel_pos;
      ci.gps.vel = gps_vel;

      // trajectory lerp
      // ABORT: through the read pointers, on the trajectory's own time base (traj_t0 = 0 for the nominal trajectory)
      float traj_t = timer / 1000000.0 - traj_t0;
      for (int i = 0; i < 15; i++) {
        ci.X_trg[i] = lerp(traj_states[pt][i], traj_states[pt + 1][i], traj_time[pt], traj_time[pt + 1], traj_t);
      }
      for (int i = 0; i < 4; i++) {
        ci.U_ff[i] = lerp(traj_inputs[pt][i], traj_inputs[pt + 1][i], traj_time[pt], traj_time[pt + 1], traj_t);
      }

      ci.GND_val = !flight_armed;

      if (ci.GND_val) { // GPS lock
        ci.gps.pos = {north : 0, west : 0, up : 0.31};
        ci.gps.vel = {north : 0, west : 0, up : 0};
        ci.new_gps_packet = true;
      } else if (gps_rel_pos.north != last_gps_pos.north || gps_rel_pos.west != last_gps_pos.west || gps_rel_pos.up != last_gps_pos.up) {
        ci.new_gps_packet = true;
        last_gps_pos = gps_rel_pos;
      } else {
        ci.new_gps_packet = false;
      }

      float time_s = timer / 1000000.0;
      float loop_dT = time_s - last_time_s;
      float ideal_dT = COMMAND_INTERVAL_US / 1000000.0;

      Controller_Output co = ControllerAndEstimator::get_controller_output(ci, ideal_dT, loop_dT, &cs);
      last_co = co; // ABORT: kept for the next tick's abort trigger
      float thrust_perc;
      float diffy_perc;
      get_prop_perc(co.thrust_N, co.roll_rad_sec_squared, &thrust_perc, &diffy_perc);
      if (flight_armed) {
        Prop::set_throttle_roll(thrust_perc, diffy_perc);
        GimbalServos::setGimbalAngle(-co.gimbal_pitch_deg, co.gimbal_yaw_deg);
      }

      if (flight_armed) // we only want to log flight data, not pre-flight
      {
        TrajectoryLogger::flash_log_sensor(time_s, pt, ci, co);

        if (timer - lastlogx_est > LOG_X_EST_INTERVAL_US) {
          TrajectoryLogger::log_x_est();
          lastlogx_est = timer;
        }
      }

      if (timer - lasttelemetry > TELEMETRY_INTERVAL_US) {
        lasttelemetry = timer;

        fp.ci = ci;
        fp.ci.imu.accel_x = cs.dnf_out[0];
        fp.ci.imu.accel_y = cs.dnf_out[1];
        fp.ci.imu.accel_z = cs.dnf_out[2];
        fp.ci.imu.gyro_yaw = cs.dnf_out[3];
        fp.ci.imu.gyro_pitch = cs.dnf_out[4];
        fp.ci.imu.gyro_roll = cs.dnf_out[5];
        fp.ci.mag.mag_x = cs.dnf_out[6];
        fp.ci.mag.mag_y = cs.dnf_out[7];
        fp.ci.mag.mag_z = cs.dnf_out[8];
        fp.x_est = cs.x_est;
        fp.co = co;

        fp.elapsed_time = timer / 1000000.0;
        fp.flight_armed = flight_armed;
        fp.thrust_perc = thrust_perc;
        fp.diffy_perc = diffy_perc;
        fp.rtk_status = (GPS::ubx.pvt_solution.data->flags >> 6) & 0b11;
        float hor, ver;
        int sats;
        GPS::get_gps_precision(&hor, &ver, &sats);
        fp.gps_hor_prec = hor;
        fp.gps_ver_prec = ver;
        fp.gps_sat_count = sats;

        CommandRouter::send_command("tr", fp);
      }
      counter++;

      unsigned long target_slp = COMMAND_INTERVAL_US - (timer - lastloop);
      delayMicroseconds(target_slp < COMMAND_INTERVAL_US ? target_slp : 0); // don't delay for too long
      lastloop += COMMAND_INTERVAL_US;
      last_time_s = time_s;
    }
  }

  digitalWrite(PIN_TEST_LED, HIGH);

  Logging::disarm();

  Prop::stop();

  fp.flight_armed = false;
  CommandRouter::send_command("tr", fp);

  CommsSerial.printf("Finished %ld loop iterations.\n", counter);
  if (kill_flag) {
    CommsSerial.println("Flight was terminated due to manual kill.");
  }
  if (auto_land) {
    CommsSerial.println("Autoland sequence was activated.");
  }
}

void heartbeat() {
  last_hb = timer;
}

// add relevant router cmds
void begin() {
  CommandRouter::add(start_flight_loop, "start_flight_loop");
  CommandRouter::add(heartbeat, "hb", "send every 0.5 seconds or flight will terminate");
  CommandRouter::add_flag(&kill_flag, "k", "terminate the flight loop early");
  CommandRouter::add_flag(&auto_land_cmd, "a", "trigger auto land sequence");
  CommandRouter::add_flag(&arm_flag, "arm", "start following a trajectory");
}

// prompt user for log file name, then follow trajectory
void start_flight_loop() {
  if (!Logging::is_armed()) {
    CommsSerial.println("FAILURE: Flash logging is not armed. Use command log_arm.");
    return;
  }

  CommsSerial.println("Flight loop starting. Type arm to start following trajectory.");
  follow_trajectory();
  CommsSerial.println("Finished following trajectory!");
}

} // namespace TrajectoryFollower
