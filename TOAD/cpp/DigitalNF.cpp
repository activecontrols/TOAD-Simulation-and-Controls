// Digital Sequential Notch Filter implementation for ASTRAv2.
// All IMU sensors are meant to be passed to this array always. Includes a
// GND switch option so the filter does not run on the ground when no
// vibrations are present.

#include "math.h"
#include "matlab_funcs.h"

Vector9 FilterStage(Vector9 IN, float f0, float fs, float width, int stage, Matrix9_4 &X, Matrix9_4 &Y) {
    float w0 = 2 * M_PI * f0 / fs;
    float r = exp(-M_PI * width / fs);

    // Calculate Normalization Gain
    float G_num = 1 - 2 * r * cos(w0) + r * r;
    float G_den = 2 - 2 * cos(w0);
    float G = G_num / G_den;

    // Filter Out
    Vector9 NUM = G * (IN - 2 * cos(w0) * X.block<9, 1>(0, 2 * stage) + X.block<9, 1>(0, 2 * stage + 1));
    Vector9 DEN = 2 * r * cos(w0) * Y.block<9, 1>(0, 2 * stage) - r * r * Y.block<9, 1>(0, 2 * stage + 1);
    Vector9 OUT = NUM + DEN;

    // Update memory
    X.block<9, 1>(0, 2 * stage + 1) = X.block<9, 1>(0, 2 * stage);
    X.block<9, 1>(0, 2 * stage) = IN;
    Y.block<9, 1>(0, 2 * stage + 1) = Y.block<9, 1>(0, 2 * stage);
    Y.block<9, 1>(0, 2 * stage) = OUT;

    return OUT;
}

Vector9 DigitalNF(Vector9 IN, bool GND, float THRUST_PERC, float dT, Matrix9_4 &X, Matrix9_4 &Y) {
    Vector9 OUT;

    Matrix2_2 TRACK;
    TRACK << 1.3525, 42.6278, 2.6867, 84.8000;

    float fs = 1 / dT;
    float width1 = 40; // Hz

    // Sequential Notch Filter
    if (!GND) {
        // Notch #1
        // Setup second notch following track #1
        float f0 = TRACK(0, 0) * THRUST_PERC + TRACK(0, 1); // Hz
        IN = FilterStage(IN, f0, fs, width1, 0, X, Y);

        // 2nd Order Butterworth Filter
        // Set up a 2nd Order Butterworth Filter
        f0 = 12;
        float C = tan(M_PI * f0 / fs);
        float A1 = 2 * (C * C - 1) / (1 + sqrt(2) * C + C * C);
        float A2 = (1 - sqrt(2) * C + C * C) / (1 + sqrt(2) * C + C * C);
        float B0 = C * C / (1 + sqrt(2) * C + C * C);
        float B1 = 2 * B0;
        float B2 = B0;
        OUT = B0 * IN + B1 * X.block<9, 1>(0, 2) + B2 * X.block<9, 1>(0, 3) - A1 * Y.block<9, 1>(0, 2) - A2 * Y.block<9, 1>(0, 3);

        X.block<9, 1>(0, 3) = X.block<9, 1>(0, 2);
        X.block<9, 1>(0, 2) = IN;
        Y.block<9, 1>(0, 3) = Y.block<9, 1>(0, 2);
        Y.block<9, 1>(0, 2) = OUT;
    } else {
        // The output is just the input(passthrough)
        OUT = IN;

        // --- Update states sequentially for a smooth switch-on ---

        // Stage 1: History is updated with IN and OUT
        X.block<9, 1>(0, 1) = X.block<9, 1>(0, 0);
        X.block<9, 1>(0, 0) = IN;
        Y.block<9, 1>(0, 1) = Y.block<9, 1>(0, 0);
        Y.block<9, 1>(0, 0) = OUT;

        // The input to stage 2 is the output of stage 1
        Vector9 IN_2 = OUT;
        Vector9 OUT_2 = IN_2; // Passthrough for stage 2

        // Stage 2: History is updated
        X.block<9, 1>(0, 3) = X.block<9, 1>(0, 2);
        X.block<9, 1>(0, 2) = IN_2;
        Y.block<9, 1>(0, 3) = Y.block<9, 1>(0, 2);
        Y.block<9, 1>(0, 2) = OUT_2;

        // The final output is the passthrough from the last stage
        OUT = OUT_2;
    }

    return OUT;
}