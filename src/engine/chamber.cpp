#include "chamber.hpp"
#include <QtConcurrent/QtConcurrent>

namespace phys {

void Chamber::fillRandom(size_t N, VelocityVal maxV, Mass m, Length r) {
    for (size_t i = 0; i < N; ++i) {
        Velocity v = randomSphere<Unit<num_t>>() * maxV;
        v *= randomShift();
        m_atoms.push_back(GasAtom{randomInCube(m_chamberCorner) *= 0.8, v, m, r});
    }
}

void Chamber::fillRandomAxis(size_t N, VelocityVal maxV, Mass m, Length r, size_t axis) {
    for (size_t i = 0; i < N; ++i) {
        Velocity v{};
        v[axis] = maxV * randomShift();
        m_atoms.push_back(GasAtom{randomInCube(m_chamberCorner) *= 0.9, v, m, r});
    }
}

void Chamber::step() {
    m_atoms.move(m_dt);

    // for (size_t i = 0; i < m_atoms.size(); ++i) {
    //     handleWallCollision(i);
    // }

    m_atoms.handleWallCollisions();

    if (m_enableCollision) {
#if 0
        for (size_t i = 0; i < m_atoms.size(); ++i) {
            for (size_t j = i + 1; j < m_atoms.size(); ++j) {
                handleCollision(i, j);
            }
        }
#else 
        m_atoms.handleCollisions([this] (size_t i, size_t j) {
            handleCollision(i, j);
        });
#endif
    }
    m_time += m_dt;
}

void Chamber::getMetrics(Metrics& metrics) const {
    metrics.chamberCorner = m_chamberCorner;

    metrics.kineticEnergy = Energy{};
    metrics.atoms.resize(m_atoms.size());

    for(size_t i = 0; i < m_atoms.size(); ++i) {
        metrics.atoms[i] = m_atoms.getAtom(i);
        metrics.kineticEnergy += metrics.atoms[i].getKineticDistributed();
    }

    metrics.time = m_time;
    metrics.volume = Volume{1.};

    for (size_t i = 0; i < UniverseDim; ++i) {
        metrics.volume *= *m_chamberCorner[i]; // HACK: I sozdal. I ignore.
    }

    for (size_t i = 0; i < 2 * UniverseDim; ++i) {
        metrics.pressure[i] = m_atoms.getWallImpulse(i) / (m_time - m_impulseMeasureStart) /
                              (metrics.volume / m_chamberCorner[i / 2]);
        if (metrics.pressure[i] < Pressure{0.})
            metrics.pressure[i] *= -1.;
    }
}

bool Chamber::hasCollision(size_t i, size_t j) {
    return (m_atoms[i].getPos() - m_atoms[j].getPos()).Len2() <
           (m_atoms[i].getRadius() + m_atoms[j].getRadius()) *
               (m_atoms[i].getRadius() + m_atoms[j].getRadius());
}

void Chamber::handleCollision(size_t i, size_t j) {
    if (!hasCollision(i, j))
        return;
    m_atoms[i].collide(m_time);
    m_atoms[j].collide(m_time);

    Velocity v1 = m_atoms[i].getVelocity();
    Velocity v2 = m_atoms[j].getVelocity();

    Mass m1 = m_atoms[i].getMass();
    Mass m2 = m_atoms[j].getMass();

    Position r1 = m_atoms[i].getPos();
    Position r2 = m_atoms[j].getPos();

    auto axis = Normalize(r2 - r1);

    VelocityVal pj1 = (v1, axis);
    VelocityVal pj2 = (v2, axis);

    if(pj1 < pj2) {
        return;
    }

    v1 -= axis * pj1;
    v2 -= axis * pj2;

    VelocityVal nex_pj1 = ((m2* pj2 *= 2.) + pj1 * (m1 - m2)) / (m1 + m2);
    VelocityVal nex_pj2 = ((m1* pj1 *= 2.) + pj2 * (m2 - m1)) / (m1 + m2);

    v1 += axis * nex_pj1;
    v2 += axis * nex_pj2;

    m_atoms[i].setVelocity(v1);
    m_atoms[j].setVelocity(v2);
}

void Chamber::handleWallCollision(size_t i) {
    Position r = m_atoms[i].getPos();
    Velocity v = m_atoms[i].getVelocity();

    for (size_t j = 0; j < UniverseDim; ++j) {
        if (r[j] < m_atoms[i].getRadius()) {
            r[j] = ((m_atoms[i].getRadius()) * num_t{2}) - r[j];
            v[j] = -v[j];
            m_wallImpulse[2 * j] += m_atoms[i].getMass() * v[j] *= 2;
            m_atoms[i].setPos(r);
            m_atoms[i].setVelocity(v);
        } else if (r[j] + m_atoms[i].getRadius() > m_chamberCorner[j]) {
            r[j] = ((m_chamberCorner[j] - m_atoms[i].getRadius()) * num_t{2}) - r[j];
            v[j] = -v[j];
            m_wallImpulse[2 * j + 1] += m_atoms[i].getMass() * v[j] *= 2;
            m_atoms[i].setPos(r);
            m_atoms[i].setVelocity(v);
        }
    }
}

} // namespace phys
