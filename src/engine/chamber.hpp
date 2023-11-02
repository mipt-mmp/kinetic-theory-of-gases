#ifndef ENGINE_UNIVERSE_HPP
#define ENGINE_UNIVERSE_HPP

#include "ballsCollection.hpp"
#include "gasAtom.hpp"

namespace phys {

const Time defaultDeltaTime = 1e4_sec;
const size_t AtomsPerCell = 10;

class Chamber {
    Position m_chamberCorner;
    // std::vector<GasAtom> m_atoms;
    BallsCollection m_atoms;
    Time m_time;
    Time m_dt = 0.01_sec;

    bool m_enableCollision = true;

    static constexpr const Time min_dt = 1e-1_sec;
    static constexpr const Time max_dt = 1e4_sec;

    Time m_impulseMeasureStart;
    std::array<phys::ImpulseVal, 6> m_wallImpulse;


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
    Chamber(Position corner = {})
        : m_chamberCorner(corner), m_atoms(std::max(corner.X(), corner.Y()), 1_sec) {
            m_atoms.setWalls(corner);
            m_atoms.setCellSize(1e-9_m);
        }

    void fillRandom(size_t N, VelocityVal maxV, Mass m, Length r);

    void fillRandomAxis(size_t N, VelocityVal maxV, Mass m, Length r, size_t axis = 0);

    void fillRandomHalf(size_t N, VelocityVal maxV, Mass m, Length r, int half);

    void updateCellSize();

    void setWalls(Position pos) {
            m_chamberCorner = pos;
            m_atoms.setWalls(pos);
    }

    void step();

    void getMetrics(Metrics& metrics) const;

    void setDT(Time dt) {
        m_dt = dt;
    }

    void setXLength(Length len);

    void openHole(bool open) {
        m_atoms.setEnableHole(open);
    }

private:
    bool hasCollision(size_t i, size_t j);

    void handleCollision(size_t i, size_t j);

    void handleWallCollision(size_t i);
};

} // namespace phys

#endif /* ENGINE_UNIVERSE_HPP */
