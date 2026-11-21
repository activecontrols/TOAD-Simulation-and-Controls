// Verification only (never part of the flight build): host solve time of Abort::plan over the verification abort starts.
// Usage: abort_timing results/abort_cases_X_est.csv   (written by verify_polyqp.m: one 15-element X_est per row)
// Each start is planned with the default command (hover thrust, gimbal centred, no disturbance), as in verify_polyqp.
// Each case is timed 5 times and the fastest run is kept (that removes scheduler noise, not cold-cache effects:
// the first run of each case is also reported as the worst first-call time).
#include "Abort.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

static float TIME[Abort::MAX_POINTS];
static float STATES[Abort::MAX_POINTS][15];
static float INPUTS[Abort::MAX_POINTS][4];

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s abort_cases_X_est.csv\n", argv[0]);
    return 2;
  }
  std::FILE *f = std::fopen(argv[1], "r");
  if (f == nullptr) {
    std::fprintf(stderr, "cannot open %s\n", argv[1]);
    return 2;
  }
  std::vector<Vector15> X;
  for (;;) {
    double v[15];
    int n = 0;
    for (; n < 15; ++n) {
      if (std::fscanf(f, n == 0 ? " %lf" : " ,%lf", &v[n]) != 1) break;
    }
    if (n != 15) break;
    Vector15 x;
    for (int i = 0; i < 15; ++i) x(i) = static_cast<float>(v[i]);
    X.push_back(x);
  }
  std::fclose(f);

  const Vector4 hover(0.0f, 0.0f, 1.1690f * 9.80145f, 0.0f);
  const Vector3 zero3 = Vector3::Zero();
  std::vector<double> best_ms, first_ms;
  int counts[3] = {0, 0, 0};
  for (const Vector15 &x : X) {
    double best = 1e30;
    Abort::Result r{};
    for (int rep = 0; rep < 5; ++rep) {
      const auto t0 = std::chrono::steady_clock::now();
      r = Abort::plan(x, hover, zero3, TIME, STATES, INPUTS, Abort::MAX_POINTS);
      const auto t1 = std::chrono::steady_clock::now();
      const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
      if (rep == 0) first_ms.push_back(ms);
      best = std::min(best, ms);
    }
    best_ms.push_back(best);
    ++counts[static_cast<int>(r.status)];
  }
  if (best_ms.empty()) {
    std::fprintf(stderr, "no cases read\n");
    return 2;
  }
  std::vector<double> sorted = best_ms;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted.size() % 2 ? sorted[sorted.size() / 2] : 0.5 * (sorted[sorted.size() / 2 - 1] + sorted[sorted.size() / 2]);
  const size_t worst_i = static_cast<size_t>(std::max_element(best_ms.begin(), best_ms.end()) - best_ms.begin());
  std::printf("Abort::plan host solve time over %zu verification starts (OK %d, INFEASIBLE_FLOWN %d, FALLBACK %d):\n", best_ms.size(), counts[0], counts[1],
              counts[2]);
  std::printf("  median %.3f ms, p95 %.3f ms, worst %.3f ms (case %zu, z0 = %.2f m); worst first call %.3f ms\n", median,
              sorted[static_cast<size_t>(0.95 * static_cast<double>(sorted.size() - 1))], best_ms[worst_i], worst_i + 1, static_cast<double>(X[worst_i](6)),
              *std::max_element(first_ms.begin(), first_ms.end()));
  return 0;
}
