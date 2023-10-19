#include "universe.hpp"
#include "physconstants.hpp"

namespace phys {

Time Universe::getTime() const
{
    return m_time;
}

void Universe::setFluctoationPerion(const Time& newFluctoationPerion)
{
    m_fluctoationPerion = newFluctoationPerion;
}

Distance Universe::getMassCenter() const {
    MassMoment mass_moments;
    Mass mass_sum;
    
    for (auto* mp : m_mps) {
        mass_moments += mp->getPos() * mp->getMass();
        mass_sum += mp->getMass();
    }

    return mass_moments / mass_sum;
}

Velocity Universe::getVelocityCenter() const {
    Impulse sum_impulse{};
    Mass mass_sum{};

    for (auto* mp : m_mps) {
        sum_impulse += mp->getVelocity() * mp->getMass();
        mass_sum    += mp->getMass();
    }

    return sum_impulse / mass_sum;
}

void Universe::fluctuate(num_t degree)
{
    auto* mp = m_mps[rand() % m_mps.size()];
    Force f = phys::random<ForceVal>();
    f *= degree;
    mp->applyForce(f);
}

void Universe::recalcOptimalDt() {

#if 0
    VelocityVal max_speed{};

    for (auto* mp : m_mps) {
        max_speed = std::max(mp->getVelocity().Len(), max_speed);
    }

    const LengthVal normalizer = 2.5e3_m;
    m_dt = normalizer / max_speed;
#endif

    num_t koef = INFINITY;

    for(const auto* mp : m_mps) {
        for(const auto* oth : m_mps) {
            if(oth == mp) continue;
            koef = std::min(koef, *(mp->getPos()-oth->getPos()).Len() / *(mp->getVelocity() - oth->getVelocity()).Len());
        }
    }

    m_dt = (1_sec *= koef * config::DT_NORMALIZER);
    m_dt = std::min(max_dt, std::max(min_dt, m_dt));
}

Time Universe::getOptimalDt() {
    return m_dt;
}

Length Universe::getMaxDist() const
{
    Length maxLen = 0_m;
    for(size_t i = 0; i < m_mps.size(); ++i) {
        maxLen = std::max(maxLen, m_mps[i]->getPos().Len());
    }
    return maxLen;
}

void Universe::shiftMassCenter()
{
    Distance new_center = getMassCenter();
    Velocity shift_velocity = getVelocityCenter();

    for (auto* mp : m_mps) {
        mp->setPos(mp->getPos() - new_center);
        mp->setVelocity(mp->getVelocity() - shift_velocity);
    }
}

void Universe::simulateStep(Time dt) {
    applyGravitation();
    
    if((m_fluctoationPerion <=> 0_sec) != std::partial_ordering::equivalent && (m_time - m_lastFluctoation) > m_fluctoationPerion) {
        fluctuate(config::FLUCTATION_DEGREE);
        m_lastFluctoation = m_time;
    }

    recalcOptimalDt();
    dt = m_dt;
    for (auto* mp: m_mps) {
        mp->move(dt);
    }

    fixUniverse(getMetrics());

    m_time += dt;
}

void Universe::fixUniverse(const Metrics& cur_metrics) {
    // Fix works only for 2 points system :D
    if (m_mps.size() != 2) {
        return;
    }

    static const Metrics begin_metrics = cur_metrics;
    
    // new_p - p = delta_p
    // 
    
    VelocityVal max_velocity{-INFINITY};
    MaterialPoint* fastest_point = nullptr;
    for (auto* mp: m_mps) {
        if (mp->getVelocity().Len() > max_velocity) {
            max_velocity = mp->getVelocity().Len();
            fastest_point = mp;
        }
    }
    assert(fastest_point);

    EnergyVal delta_energy = cur_metrics.energy - begin_metrics.energy;
    // a^2mv^2/2 = mv^2/2 - c
    // a^2 = 1 + (2c)/(mv^2)
    Unit<num_t> coef = Unit<num_t>{1.0} - (Unit<num_t>{2.0} * delta_energy) / (fastest_point->getMass() * max_velocity * max_velocity);

    Velocity new_velocity = fastest_point->getVelocity();
    new_velocity *= (*coef);
    fastest_point->setVelocity(new_velocity);

    auto get_errors = [] (const Metrics& first, const Metrics& second) {
        return std::make_pair((first.impulse - second.impulse).Len(),
                              (first.impulsemoment - second.impulsemoment).Len());
    };

    std::pair<ImpulseVal, ImpulseMomentVal> best_errors{get_errors(getMetrics(), begin_metrics)};

    Velocity unrotated_velocity = fastest_point->getVelocity();
    double degree_step = M_PI / 45.0;
    for (double cur_degree = degree_step; cur_degree < 2.0 * M_PI; cur_degree += degree_step) {
        auto cosine = Unit<num_t>{cos(cur_degree)};
        auto sinus = Unit<num_t>{sin(cur_degree)};

        Velocity rotated_velocity{unrotated_velocity.m_coord[0] * cosine - unrotated_velocity.m_coord[1] * sinus,
                                  unrotated_velocity.m_coord[0] * sinus + unrotated_velocity.m_coord[1] * cosine};

        fastest_point->setVelocity(rotated_velocity);
        auto rotated_errors = get_errors(getMetrics(), begin_metrics);

        if ((rotated_errors.first < best_errors.first) &&
            (rotated_errors.second < best_errors.second)) {
            best_errors = rotated_errors;
            new_velocity = rotated_velocity;
        }
    }

    fastest_point->setVelocity(new_velocity);
}

void Universe::applyGravitation() {
    for(size_t i = 0; i < m_mps.size(); ++i) {
        for(size_t j = 0; j < m_mps.size(); ++j) {
            if(i == j) continue;
            auto direction = -Normalize(m_mps[i]->getPos() - m_mps[j]->getPos());
            assert(std::abs(*direction.Len2() - 1) < 1e-4l);
            auto dist      =          (m_mps[i]->getPos() - m_mps[j]->getPos()).Len2();
            assert(*dist > 0.l);
            Force f = direction * phys::consts::G * m_mps[i]->getMass() * m_mps[j]->getMass() / dist;
            m_mps[i]->applyForce(f);
        }
    }
}

ImpulseMoment Universe::getImpulseMoment() const
{
    ImpulseMoment p{};
    for(const auto* mp: m_mps) {
        p += mp->getImpulseMoment();
    }

    return p;
}

Impulse Universe::getImpulse() const
{
    Impulse p{};
    for(const auto* mp: m_mps) {
        p += mp->getImpulse();
    }

    return p;
}

Energy Universe::getEnergy() const {
    Energy e{};
    for(const auto* mp: m_mps) {
        e += mp->getKinetic();
    }
    e += getPotentialGravitationEnergy();

    return e;
}

Energy Universe::getPotentialGravitationEnergy() const {
    Energy e;
    for(size_t i = 0; i < m_mps.size(); ++i) {
        for(size_t j = i+1; j < m_mps.size(); ++j) {
            auto dist = (m_mps[i]->getPos() - m_mps[j]->getPos()).Len();
            e += -phys::consts::G * m_mps[i]->getMass() * m_mps[j]->getMass() / dist;
        }
    }
    return e;
}



}
