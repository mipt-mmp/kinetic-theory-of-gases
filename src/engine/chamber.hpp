#ifndef ENGINE_UNIVERSE_HPP
#define ENGINE_UNIVERSE_HPP

#include "gasAtom.hpp"
namespace phys {

const Time defaultDeltaTime = 1e4_sec;

class Chamber {
  std::vector<GasAtom> m_atoms;
  Time m_time;
  Time m_dt = 0.01_sec;

  bool m_enableCollision = true;

  static constexpr const Time min_dt = 1e-1_sec;
  static constexpr const Time max_dt = 1e4_sec;
  
  Time m_impulseMeasureStart;
  std::array<phys::ImpulseVal, 6> m_wallImpulse;
  
  Position m_chamberCorner;

public:

  enum ChamberWall {
    Left = 0,
    Right,
    Top,
    Bottom,
    Front,
    Back,
  };

  struct Metrics {
      Position chamberCorner;
      std::vector<GasAtom> atoms;
      Volume volume;
      Vector<Energy> kineticEnergy;
      std::array<Pressure, 2 * UniverseDim> pressure;
      Impulse impulse;
      ImpulseMoment impulseMoment;
      Time time;
  };
public:
  
  Chamber(Position corner) : m_chamberCorner(corner) {}

  void fillRandom(size_t N, VelocityVal maxV, Mass m, Length r);

  void fillRandomAxis(size_t N, VelocityVal maxV, Mass m, Length r, size_t axis = 0);

  void step();

  void getMetrics(Metrics& metrics) const;

  void setDT(Time dt) {m_dt = dt;}

private:

  bool hasCollision(size_t i, size_t j) const;

  void handleCollision(size_t i, size_t j);

  void handleWallCollision(size_t i);


};

}

#endif /* ENGINE_UNIVERSE_HPP */
