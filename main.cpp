#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include "cpgplot.h"

struct Wheel {
    int nbuckets;
    double damping;
    double inertia;
    double drainrate;
    double fillrate;
    double gravity;
    double radius;
    double rotation;
    double velocity;
    std::vector<double> buckets;
    std::vector<std::pair<float, float>> com;
};

struct Derivative {
    double rdot;
    double vdot;
    std::vector<double> bdot;
};

Wheel createWheel(int n) {
    return Wheel {
        .nbuckets = n,
        .damping = 3.0,
        .inertia = 0.1,
        .drainrate = 0.3,
        .fillrate = 0.3,
        .gravity = 9.5,
        .radius = 1.0,
        .rotation = 0.5 * 2 * M_PI,
        .velocity = 0.4,
        .buckets = std::vector<double>(n, 0.0),
        .com = std::vector<std::pair<float, float>>()
    };
}

Derivative derive(const Wheel& wheel) {
    const int nbuckets = wheel.nbuckets;
    const double radius = wheel.radius;
    const double gravity = wheel.gravity;
    const double inertia = wheel.inertia + wheel.radius * wheel.radius * std::accumulate(wheel.buckets.begin(), wheel.buckets.end(), 0.0);

    const double dampingVelocity = -wheel.damping * wheel.velocity;
    const double rg = radius * gravity;
    std::vector<double> bdot(nbuckets, 0.0);
    const double f = wheel.fillrate / 2;

    const double rotationIncrement = 2 * M_PI / nbuckets;
    const double cos_rotationIncrement = cos(rotationIncrement);

    double torque = dampingVelocity;

    for (int i = 0; i < nbuckets; ++i) {
        const double r = wheel.rotation + i * rotationIncrement;
        const double sin_r = sin(r);
        const double cos_r = cos(r);

        torque += rg * wheel.buckets[i] * sin_r;

        bdot[i] = -wheel.drainrate * wheel.buckets[i];

        if (cos_r > abs(cos_rotationIncrement)) {
            const double x = atan2(tan(r), 1);
            bdot[i] += f * (cos(nbuckets * x / 2) + 1);
        }
    }

    return {
        .rdot = wheel.velocity,
        .vdot = (nbuckets > 0) ? torque / inertia : 0.0,
        .bdot = bdot
    };
}

Wheel clone(const Wheel& wheel) {
    Wheel dup = wheel;
    dup.buckets = wheel.buckets;
    return dup;
}

void apply(Wheel& wheel, const Derivative& wdot, double dt) {
    wheel.rotation += wdot.rdot * dt;
    wheel.velocity += wdot.vdot * dt;
    for (int i = 0; i < wheel.nbuckets; ++i) {
        wheel.buckets[i] += wdot.bdot[i] * dt;
    }
}

Derivative rk4(const Wheel& k1, double dt) {
    const Derivative k1d = derive(k1);
    Wheel k2 = clone(k1);
    apply(k2, k1d, dt * 0.5);
    const Derivative k2d = derive(k2);
    Wheel k3 = clone(k1);
    apply(k3, k2d, dt * 0.5);
    const Derivative k3d = derive(k3);
    Wheel k4 = clone(k1);
    apply(k4, k3d, dt);
    const Derivative k4d = derive(k4);

    Derivative dot = k1d;
    dot.rdot += 2 * (k2d.rdot + k3d.rdot) + k4d.rdot;
    dot.vdot += 2 * (k2d.vdot + k3d.vdot) + k4d.vdot;

    for (int i = 0; i < k1.nbuckets; ++i) {
        dot.bdot[i] += 2 * (k2d.bdot[i] + k3d.bdot[i]) + k4d.bdot[i];
    }

    Wheel result = clone(k1);
    apply(result, dot, dt / 6);
    return derive(result);
}

std::pair<float, float> calculate_com(const Wheel& wheel) {
    double total_mass = 0.0;
    double x_sum = 0.0;
    double y_sum = 0.0;

    for (int i = 0; i < wheel.nbuckets; ++i) {
        total_mass += wheel.buckets[i];

        double r = wheel.rotation + i * 2 * M_PI / wheel.nbuckets;
        x_sum += sin(r) * wheel.buckets[i];
        y_sum -= cos(r) * wheel.buckets[i];
    }

    float x_cm = 0.0f;
    float y_cm = 0.0f;
    if (total_mass != 0.0) {
        x_cm = static_cast<float>(x_sum / total_mass);
        y_cm = static_cast<float>(y_sum / total_mass);
    }

    return std::make_pair(x_cm, y_cm);
}

void simulateAndPlot() {
    if (cpgbeg(0, "/xwindow", 1, 1) != 1) {return;}
    
    cpgenv(-1, 1, -1, 1, 0, 0);
    cpglab("x", "y", "Center of Mass of Water Wheel");

    Wheel wheel = createWheel(15);
    double dt = 0.01;
    double t = 0.0;

    std::vector<float> x_cm_data;
    std::vector<float> y_cm_data;

    while (true) { 
        // replace with while(true) for infinite loop
        // replace with t<n for finite loop
        std::pair<float, float> cm = calculate_com(wheel);
        x_cm_data.push_back(cm.first);
        y_cm_data.push_back(cm.second);

        cpgsci(2);
        cpgpt(x_cm_data.size(), &x_cm_data[0], &y_cm_data[0], 1);

        Derivative wdot = rk4(wheel, dt);
        apply(wheel, wdot, dt);

        t += dt;
    }
    printf("Done\n");
    cpgend();
}

int main() {
    simulateAndPlot();
    return 0;
}
