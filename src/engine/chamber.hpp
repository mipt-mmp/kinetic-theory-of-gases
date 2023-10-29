#ifndef ENGINE_UNIVERSE_HPP
#define ENGINE_UNIVERSE_HPP

#include "gasAtom.hpp"
#include "sync.hpp"

#include <QAtomicInteger>
#include <QThreadPool>

namespace phys {

class Chamber;

struct ChamberBlock {
    size_t blockStart;
    size_t blockEnd;
};

class ChamberBlockRunner : public QRunnable {
public:
    explicit ChamberBlockRunner(Chamber& chamber);
    ChamberBlockRunner(Chamber& chamber, uint64_t blockId);

    uint64_t getBlockId() const;
    void setBlockId(uint64_t blockId);

    void run() final;

private:
    void handleCollisionsForBlock(size_t blockId);

private:
    uint64_t m_blockId;
    Chamber& m_chamber;
};

class ChamberMutexUnlocker : public QRunnable {
public:
    explicit ChamberMutexUnlocker(Chamber& chamber);

    void run() final;

private:
    Chamber& m_chamber;
};

class Chamber {
    friend class ChamberBlockRunner;
    friend class ChamberMutexUnlocker;

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
    explicit Chamber(Position corner);

    void fillRandom(size_t N, VelocityVal maxV, Mass m, Length r);
    void fillRandomAxis(size_t N, VelocityVal maxV, Mass m, Length r, size_t axis = 0);

    void step();

    void getMetrics(Metrics& metrics) const;
    size_t getBlocksNumber() const;

    void setDt(Time dt);

private:
    bool hasCollision(size_t i, size_t j) const;
    void handleCollision(size_t i, size_t j);
    void handleWallCollision(size_t i);

    void updateBlockForAllAtoms();
    void updateBlockForSingleAtom(size_t atomId);
    void updateBlockSegments();

private:
    static const size_t k_maxDimensions = 3;
    static const size_t k_maxThreads = 8;
    static const size_t k_blocksPerDim = 64;

    std::vector<GasAtom> m_atoms;

    Time m_time;
    Time m_dt{0.01_sec};

    bool m_enableCollisions{true};

    static constexpr const Time min_dt{1e-1_sec};
    static constexpr const Time max_dt{1e4_sec};

    Time m_impulseMeasureStart;
    std::array<phys::ImpulseVal, k_maxDimensions * 2> m_wallImpulse;

    Position m_chamberCorner;

    std::vector<ChamberBlock> m_blocks;

    QMutex m_mutex;
    QSync::QBarrier m_blockBarrier;
    QThreadPool m_threadPool;
};

} // namespace phys

#endif /* ENGINE_UNIVERSE_HPP */
