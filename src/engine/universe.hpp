#ifndef ENGINE_UNIVERSE_HPP
#define ENGINE_UNIVERSE_HPP

#include "materialpoint.hpp"
namespace phys {

const Time defaultDeltaTime = 1e4_sec;

class Universe {
  std::vector<MaterialPoint *> m_mps;
  Time m_time;
  Time m_dt;

  Time m_lastFluctoation = 0_sec;
  Time m_fluctoationPerion = 0_sec;


  static constexpr const Time min_dt = 1e-1_sec;
  static constexpr const Time max_dt = 1e4_sec;

public:
  struct Metrics {
      Impulse impulse;
      Energy energy;
      ImpulseMoment impulsemoment;
      Time time;
  };

  Metrics getMetrics() const {
      return {
        getImpulse(),
        getEnergy(),
        getImpulseMoment(),
        m_time
      };
  }

  void addMaterialPoint(MaterialPoint *mp) { m_mps.push_back(mp); }

  Position getMassCenter() const;
  Velocity getVelocityCenter() const;

  void fluctuate(num_t degree);

  void recalcOptimalDt();
  Time getOptimalDt();

  Length getMaxDist() const;

  void shiftMassCenter();

  void simulateStep(Time dt = defaultDeltaTime);

  void fixUniverse(const Metrics& cur_metrics);

  ImpulseMoment getImpulseMoment() const;
  Impulse getImpulse() const;
  Energy getEnergy() const;

  [[nodiscard]] Time getTime() const;

  void setFluctoationPerion(const Time& newFluctoationPerion);

  private:
  Energy getPotentialGravitationEnergy() const;

  void applyGravitation();
};

}

#endif /* ENGINE_UNIVERSE_HPP */
