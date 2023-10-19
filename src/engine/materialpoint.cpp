#include "materialpoint.hpp"

namespace phys {

void MaterialPoint::move(Time dt) {
    Acceleration a = cur_f_ / m_;
    pos_ += v_ * dt + (a * dt * dt /= 2);
    v_ += a * dt;

    prev_a  = a;
    prev_dt = dt;

    reset();
}

void MaterialPoint::applyForce(Force f) {
    cur_f_ += f;
}

void MaterialPoint::reset() {
    cur_f_ = {};
}

}
