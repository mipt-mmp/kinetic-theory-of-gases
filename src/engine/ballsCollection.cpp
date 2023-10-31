#include "ballsCollection.hpp"
#include <numeric>
#include <iostream>
#include <QtConcurrent/QtConcurrent>
#include <QMutexLocker>

static const auto StepSize = 16;
namespace phys {

static const num_t holeSize = 0.1;

detail::GasAtomProxy::GasAtomProxy(BallsCollection& balls, size_t index) : m_balls(balls), m_index(index), m_atom(balls.getAtom(index)) {}

void detail::GasAtomProxy::updateBalls() {
    for(size_t i = 0; i < UniverseDim; ++i) {
        m_balls.m_coords    [i][m_index] = *(m_atom.getPos()     [i] / m_balls.m_mScale);
        m_balls.m_velocities[i][m_index] = *(m_atom.getVelocity()[i] / m_balls.m_mScale * m_balls.m_tScale);
    }
}

void BallsCollection::deleteAtom(size_t i) {
    --m_nAtoms;

    for(size_t j = 0; j < UniverseDim; ++j) {
        std::swap(m_coords[j][i], m_coords[j][m_nAtoms]);
        m_coords[j].pop_back();

        std::swap(m_velocities[j][i], m_velocities[j][m_nAtoms]);
        m_velocities[j].pop_back();
    }

    std::swap(m_masses[i], m_masses[m_nAtoms]);
    m_masses.pop_back();

    std::swap(m_radiuses[i], m_radiuses[m_nAtoms]);
    m_radiuses.pop_back();
}

GasAtom BallsCollection::getAtom(size_t i) const {
    assert(i < m_nAtoms);
    Position pos;
    Velocity v;
    for(size_t j = 0; j < UniverseDim; ++j) {
        pos[j] = m_mScale            * m_coords    [j][i];
        v  [j] = m_mScale / m_tScale * m_velocities[j][i];
    }
    return GasAtom{pos, v, Mass{m_masses[i]}, m_mScale * m_radiuses[i]};
}

void BallsCollection::move(Time dt) {
    num_t time = *(dt / m_tScale);

    QFutureSynchronizer<void> synchronizer = {};
    for(size_t d = 0; d < UniverseDim; ++d) {
        for(size_t start = 0; start < m_nAtoms; start += m_nAtoms / StepSize) {
            synchronizer.addFuture(QtConcurrent::run( 
                [this, d, time, start] () {
                    for(size_t i = start; i < std::min(m_nAtoms, start + m_nAtoms / StepSize); ++i) {
                        m_coords[d][i] += m_velocities[d][i] * time;
                    }
                }
            ));
        }
    }
    synchronizer.waitForFinished();
}

void BallsCollection::handleWallCollisions() {
    QMutex deleteProtector = {};
    std::vector<size_t> deleteCandidates = {};

    auto isInHole = [this] (size_t atomIdx) {
        bool flag = true;
        
        for (size_t holeDim = 1; holeDim < UniverseDim; holeDim++) {
            if ((std::abs(m_coords[holeDim][atomIdx] - (m_walls[holeDim] / 2)) / m_walls[holeDim]) > holeSize) {
                flag = false;
            }
        }

        return flag;
    };
    
    QFutureSynchronizer<void> synchronizer = {};

    for (size_t j = 0; j < UniverseDim; ++j) {
        for(size_t start = 0; start < m_nAtoms; start += m_nAtoms / StepSize) {
            synchronizer.addFuture(QtConcurrent::run(
                [this, start, j, isInHole, &deleteProtector, &deleteCandidates] () {
                    for(size_t i = start; i < std::min(m_nAtoms, start + m_nAtoms / StepSize); ++i) {
                        if (m_coords[j][i] < m_radiuses[i]) {
                            if ((j == 0) && isInHole(i)) {
                                std::lock_guard guard{deleteProtector};
                                deleteCandidates.push_back(i);
                            }

                            m_coords[j][i] = (m_radiuses[i] * 2) - m_coords[j][i];
                            m_velocities[j][i] = -m_velocities[j][i];
                            m_wallImpulse[2 * j] += m_masses[i] * m_velocities[j][i] * 2;
                        } else if (m_coords[j][i] + m_radiuses[i] > m_walls[j]) {
                            m_coords[j][i] = ((m_walls[j] - m_radiuses[i]) * 2) - m_coords[j][i];
                            m_velocities[j][i] = -m_velocities[j][i];
                            m_wallImpulse[2 * j + 1] += m_masses[i] * m_velocities[j][i] * 2;
                        }
                    }
                }
            ));
        }
    }

    synchronizer.waitForFinished();

    for (auto& candidateIdx : deleteCandidates) {
        deleteAtom(candidateIdx);
    }
}

static uint32_t getShift(uint32_t x) {
    if(!x) return 1;
    x--;
    unsigned shift = 1;
    while (x >>= 1) {
        shift++;
    }

    return shift;
} 

void BallsCollection::setCellSize(Length l) {
    num_t len = *(l / m_mScale);
    m_cellSize = len;

    m_shifts[0] = 0;
    std::cout << "Grid dim: ";
    for(size_t i = 1; i < UniverseDim; ++i) {
        uint32_t nCells = std::ceil((m_walls[i-1] / len).getVal());
        std::cout << nCells << '*';
        m_shifts[i] =  m_shifts[i-1] + getShift(nCells);
    }
    uint32_t nCells = std::ceil((m_walls[UniverseDim-1] / len).getVal());
    std::cout << nCells << '\n';
    assert(m_shifts.back() + getShift(nCells) < 32);
}


void BallsCollection::handleCollisions() {
// Counting hashes;
    QFutureSynchronizer<void> syncher = {};
    for (size_t start = 0; start < m_nAtoms; start += m_nAtoms / StepSize) {
        syncher.addFuture(QtConcurrent::run(
            [start, this] () {
                for (size_t i = start; i < std::min(m_nAtoms, start + m_nAtoms / StepSize); i++) {
                    m_indicies[i] = static_cast<uint32_t>(i);
                }
            }
        ));

        syncher.addFuture(QtConcurrent::run(
            [start, this] () {
                for (size_t i = start; i < std::min(m_nAtoms, start + m_nAtoms / StepSize); i++) {
                    m_hashes[i] = 0;
                    for(size_t j = 0; j < UniverseDim; ++j) {
                        m_hashes[i] |= static_cast<uint32_t>(m_coords[j][i] / m_cellSize) << m_shifts[j];
                    }
                }
            }
        ));
    }
    syncher.waitForFinished();

    radixSort();

    m_radixBuffer.assign(m_nAtoms, 0);
    
    m_collisionList.clear();

    QFutureSynchronizer<void> synchronizer = {};

    for(size_t i = 0; i < m_nAtoms; i += m_nAtoms / StepSize) {
        synchronizer.addFuture(QtConcurrent::run(&BallsCollection::handleSub, this, i, std::min(i + m_nAtoms / StepSize, m_nAtoms)));
    }
    /*
    size_t prev = 0;
    for(size_t i = 1; i < m_nAtoms; ++i) {
        if ((m_hashes[i] != m_hashes[prev]) && ((i - prev) > 1)) {
            prev = i;
        }
    }
    synchronizer.addFuture(QtConcurrent::run(&BallsCollection::handleBlock, this, prev, m_nAtoms, action));
    */
    synchronizer.waitForFinished();
}

void BallsCollection::handleSub(size_t l, size_t r) {
    size_t prev = l;
    for(size_t i = l+1; i < r; ++i) {
        if (m_hashes[i] != m_hashes[prev]) {
            if(i > prev + 1) {
                handleBlock(prev, i);
            }
            prev = i;
        }
    }
}


void BallsCollection::radixSort() {
    const uint32_t MASK = 0xff;
    std::array<size_t, 257> sums;
    for(uint32_t shift = 0; shift < 32; shift += 8) {
        sums.fill(0);

        for(size_t i = 0; i < m_nAtoms; ++i) {
            sums[((m_hashes[i] >> shift) & MASK) + 1]++;
        }

        for(size_t i = 0; i < 256; ++i) {
            sums[i+1] += sums[i];
        }

        for(size_t i = 0; i < m_nAtoms; ++i) {
            uint32_t idx = ((m_hashes[i] >> shift) & MASK);
            m_radixIndiciesBuffer[sums[idx]  ] = m_indicies[i];
            m_radixBuffer        [sums[idx]++] = m_hashes[i];
        }

        m_radixBuffer.swap(m_hashes);
        m_radixIndiciesBuffer.swap(m_indicies);
    }
}

void BallsCollection::handleBlock(size_t l, size_t r) {
    if(r == l+2) {
        size_t i = m_indicies[l];
        size_t j = m_indicies[l+1];
        num_t dst = 0;
        for(size_t d = 0; d < UniverseDim; ++d) {
            dst += (m_coords[d][i] - m_coords[d][j]) * (m_coords[d][i] - m_coords[d][j]);
        }
        if(dst < (m_radiuses[i] + m_radiuses[j]) * (m_radiuses[i] + m_radiuses[j])) {
            QMutexLocker<QMutex> locker(&m_listMutex);
            m_collisionList.push_back(std::make_pair(i, j));
        }
        return;
    }

    for(size_t idx = l; idx < r; ++idx) {
        for(size_t jdx = idx + 1; jdx < r; ++jdx) {
            size_t i = m_indicies[idx];
            size_t j = m_indicies[jdx];
            num_t dst = 0;
            for(size_t d = 0; d < UniverseDim; ++d) {
                dst += (m_coords[d][i] - m_coords[d][j]) * (m_coords[d][i] - m_coords[d][j]);
            }
            if(dst < (m_radiuses[i] + m_radiuses[j]) * (m_radiuses[i] + m_radiuses[j])) {
                QMutexLocker<QMutex> locker(&m_listMutex);
                m_collisionList.push_back(std::make_pair(i, j));
            }
        }
    }
}


}
