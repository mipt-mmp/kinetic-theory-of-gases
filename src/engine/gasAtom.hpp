#ifndef ENGINE_MATERIALPOINT_HPP
#define ENGINE_MATERIALPOINT_HPP

#include "geometry.hpp"
#include "units.hpp"

namespace phys {

class GasAtom {

    // State characteristics
    Position m_pos;
    Velocity m_v;
    Mass m_mass;
    Length m_radius;

    decltype(Energy{} * Time{}) m_sumE;
    Time m_lastCollide;
public:
    GasAtom(const Position& pos, const Velocity& v, const Mass& m, const Length radius)
        : m_pos(pos)
        , m_v(v)
        , m_mass(m)
        , m_radius(radius) {
        if (*m_mass == 0) {
            std::cerr << "Warning: too small mass: " << m_mass << "\n";
        }
    }

    const Mass& getMass() const {
        return m_mass;
    }

    Length getRadius() const {
        return m_radius;
    }

    void setMass(const Mass& newMass) {
        m_mass = newMass;
    }

    void collide(Time t) {
        m_lastCollide = t;
    }

    Length getFreeFlight(Time t) {
        return (m_v * (t - m_lastCollide)).Len();
    }

    Energy getAverageEnergy(Time t) {
        return m_sumE / t;
    }

    const Velocity& getVelocity() const {
        return m_v;
    }

    void setVelocity(const Velocity& newVelocity) {
        m_v = newVelocity;
    }

    const Position& getPos() const {
        return m_pos;
    }

    void setPos(const Position& newPos) {
        m_pos = newPos;
    }

    void move(Time dt) {
        m_pos  += m_v * dt;
        m_sumE += getKinetic() * dt;
    }

    ImpulseMoment getImpulseMoment() const {
        auto impulse = m_v * m_mass;
        return CrossProd(m_pos, impulse);
    }

    Impulse getImpulse() const {
        return m_v * m_mass;
    }

    Energy getKinetic() const {
        Energy e = m_v.Len2() * m_mass / 2.l;
        return e;
    }

    Vector<Energy> getKineticDistributed() const {
        return MulByElement(m_v, m_v) * m_mass /= 2.l;
    }

private:
};

} // namespace phys

#endif /* ENGINE_MATERIALPOINT_HPP */
