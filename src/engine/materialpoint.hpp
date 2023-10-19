#ifndef ENGINE_MATERIALPOINT_HPP
#define ENGINE_MATERIALPOINT_HPP

#include "units.hpp"
#include "geometry.hpp"

namespace phys {

class MaterialPoint {

// State characteristics
    Position pos_;
    Velocity v_;
    Mass m_;

// Calculation only;
    Force cur_f_;
    Acceleration prev_a;
    Time prev_dt{INFINITY};
public:

    MaterialPoint(const Position& pos, const Velocity& v, const Mass& m) : pos_(pos), v_(v), m_(m) {
        if(*m_ == 0) {
            std::cerr << "Warning: too small mass: " << m_ << "\n";
        }
    }

    const Mass& getMass() const {
        return m_;
    }

    void setMass(const Mass& newMass) {
        reset();
        m_ = newMass;
    }

    const Velocity& getVelocity() const {
        return v_;
    }

    void setVelocity(const Velocity& newVelocity) {
        reset();
        v_ = newVelocity;
    }

    const Position& getPos() const {
        return pos_;
    }

    void setPos(const Position& newPos) {
        reset();
        pos_ = newPos;
    }

    void move(Time dt);

    void applyForce(Force f);

    ImpulseMoment getImpulseMoment() const {
        auto impulse = v_ * m_;

        return CrossProd(pos_, impulse);
    }

    Impulse getImpulse() const {
        return v_ * m_;
    }

    Energy getKinetic() const { 
        Energy e = v_.Len2() * m_ / 2.l; 
        return e;
    }

private:
    void reset(); 
};

}

#endif /* ENGINE_MATERIALPOINT_HPP */
