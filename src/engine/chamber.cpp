#include "chamber.hpp"

namespace phys {

ChamberBlockRunner::ChamberBlockRunner(Chamber& chamber)
    : m_chamber(chamber) {}

ChamberBlockRunner::ChamberBlockRunner(Chamber& chamber, uint64_t blockId)
    : m_chamber(chamber)
    , m_blockId(blockId) {}

uint64_t ChamberBlockRunner::getBlockId() const {
    return m_blockId;
}

void ChamberBlockRunner::setBlockId(uint64_t blockId) {
    m_blockId = blockId;
}

void ChamberBlockRunner::run() {
    auto start = m_blockId;
    auto end = m_chamber.getBlocksNumber();
    auto step = Chamber::k_maxThreads;

    for (auto i = start; i < end; i += step) {
        handleCollisionsForBlock(i);
    }

    m_chamber.m_blockBarrier.wait();
}

void ChamberBlockRunner::handleCollisionsForBlock(size_t blockId) {
    auto& block = m_chamber.m_blocks[blockId];
    bool doCollide = m_chamber.m_enableCollisions;

    auto start = block.blockStart;
    auto end = block.blockEnd;

    for (size_t i = start; i < end; i++) {
        m_chamber.handleWallCollision(i);
    }

    if (!doCollide) {
        return;
    }

    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < end; j++) {
            m_chamber.handleCollision(i, j);
        }
    }
}

ChamberMutexUnlocker::ChamberMutexUnlocker(Chamber& chamber)
    : m_chamber(chamber) {}

void ChamberMutexUnlocker::run() {
    m_chamber.m_blockBarrier.wait();
    m_chamber.m_mutex.unlock();
}

Chamber::Chamber(Position corner)
    : m_chamberCorner(corner)
    , m_blocks(std::vector<ChamberBlock>(std::pow(k_blocksPerDim, k_maxDimensions)))
    , m_blockBarrier(QSync::QBarrier(k_maxThreads + 1)) {

    // Extra 1 for ChamberMutexUnlocker
    m_threadPool.setMaxThreadCount(k_maxThreads + 1);
}

const size_t Chamber::k_maxThreads = std::thread::hardware_concurrency();

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
    m_mutex.lock();

    for (size_t i = 0; i < m_atoms.size(); ++i) {
        m_atoms[i].move(m_dt);
    }

    updateBlockForAllAtoms();
    updateBlockSegments();

    for (size_t i = 0; i < k_maxThreads; i++) {
        m_threadPool.start(new ChamberBlockRunner(*this, i));
    }

    m_threadPool.start(new ChamberMutexUnlocker(*this));
    m_time += m_dt;
}

void Chamber::getMetrics(Metrics& metrics) const {
    metrics.chamberCorner = m_chamberCorner;
    metrics.atoms = m_atoms;

    metrics.kineticEnergy = Energy{};

    for (const auto& atom : m_atoms) {
        metrics.kineticEnergy += atom.getKineticDistributed();
    }

    metrics.volume = Volume{1.};

    for (size_t i = 0; i < UniverseDim; ++i) {
        metrics.volume *= *m_chamberCorner[i]; // HACK: I sozdal. I ignore.
    }

    for (size_t i = 0; i < 2 * UniverseDim; ++i) {
        metrics.pressure[i] = m_wallImpulse[i] / (m_time - m_impulseMeasureStart) /
                              (metrics.volume / m_chamberCorner[i / 2]);

        if (metrics.pressure[i] < Pressure{0.}) {
            metrics.pressure[i] *= -1.;
        }
    }
}

size_t Chamber::getBlocksNumber() const {
    return m_blocks.size();
}

void Chamber::setDt(Time dt) {
    m_dt = dt;
}

bool Chamber::hasCollision(size_t i, size_t j) const {
    return (m_atoms[i].getPos() - m_atoms[j].getPos()).Len2() <
           (m_atoms[i].getRadius() + m_atoms[j].getRadius()) *
               (m_atoms[i].getRadius() + m_atoms[j].getRadius());
}

void Chamber::handleCollision(size_t i, size_t j) {
    if (!hasCollision(i, j)) {
        return;
    }

    Velocity v1 = m_atoms[i].getVelocity();
    Velocity v2 = m_atoms[j].getVelocity();

    Mass m1 = m_atoms[i].getMass();
    Mass m2 = m_atoms[j].getMass();

    Position r1 = m_atoms[i].getPos();
    Position r2 = m_atoms[j].getPos();

    auto axis = Normalize(r2 - r1);

    VelocityVal pj1 = (v1, axis);
    VelocityVal pj2 = (v2, axis);

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

void Chamber::updateBlockForAllAtoms() {
    for (size_t i = 0; i < m_atoms.size(); i++) {
        updateBlockForSingleAtom(i);
    }
}

void Chamber::updateBlockForSingleAtom(size_t atomId) {
    auto& atom = m_atoms[atomId];

    auto xCellSize = m_chamberCorner.X() / static_cast<double>(k_blocksPerDim);
    auto yCellSize = m_chamberCorner.Y() / static_cast<double>(k_blocksPerDim);
    auto zCellSize = m_chamberCorner.Z() / static_cast<double>(k_blocksPerDim);

    const auto k_xShift = 0;
    const auto k_yShift = 8;
    const auto k_zShift = 16;

    auto newBlock = (static_cast<uint64_t>((atom.getPos().X() / xCellSize).inner()) << k_xShift) |
                    (static_cast<uint64_t>((atom.getPos().Y() / yCellSize).inner()) << k_yShift) |
                    (static_cast<uint64_t>((atom.getPos().Z() / zCellSize).inner()) << k_zShift);

    atom.setBlockId(newBlock);
}

void Chamber::updateBlockSegments() {
    assert(!m_atoms.empty());

    std::sort(m_atoms.begin(), m_atoms.end());

    auto blockStart = 0;
    auto blockId = 0;

    auto lastId = m_atoms[0].getBlockId();

    for (size_t i = 0; i < m_atoms.size(); i++) {
        auto currId = m_atoms[i].getBlockId();

        if (currId != lastId) {
            lastId = currId;
            m_blocks[blockId].blockStart = blockStart;
            m_blocks[blockId].blockEnd = i;
            blockId++;
            blockStart = i;
        }
    }
}

} // namespace phys
