#include "ballsCollection.hpp"
#include <numeric>
#include <iostream>

namespace phys {

detail::GasAtomProxy::GasAtomProxy(BallsCollection& balls, size_t index) : m_balls(balls), m_index(index), m_atom(balls.getAtom(index)) {}

void detail::GasAtomProxy::updateBalls() {
    for(size_t i = 0; i < UniverseDim; ++i) {
        m_balls.m_coords    [i][m_index] = *(m_atom.getPos()     [i] / m_balls.m_mScale);
        m_balls.m_velocities[i][m_index] = *(m_atom.getVelocity()[i] / m_balls.m_mScale * m_balls.m_tScale);
    }
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
    for(size_t d = 0; d < UniverseDim; ++d) {
        for(size_t i = 0; i < m_nAtoms; ++i) {
            m_coords[d][i] += m_velocities[d][i] * time;
        }
    }
}

void BallsCollection::handleWallCollisions() {
    for (size_t j = 0; j < UniverseDim; ++j) {
        for(size_t i = 0; i < m_nAtoms; ++i) {
            if (m_coords[j][i] < m_radiuses[i]) {
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


void BallsCollection::countCollisions() {
// Counting hashes;
    for(size_t i = 0; i < m_nAtoms; ++i) {
        m_hashes[i] = 0;
        for(size_t j = 0; j < UniverseDim; ++j) {
            m_hashes[i] |= static_cast<uint32_t>(m_coords[j][i] / m_cellSize) << m_shifts[j];
        }
    }
    std::iota(m_indicies.begin(), m_indicies.end(), 0);

    radixSort();

    m_radixBuffer.assign(m_nAtoms, 0);
    
    m_collisionList.clear();

    size_t prev = 0;
    for(size_t i = 1; i < m_nAtoms; ++i) {
        if(m_hashes[i] != m_hashes[prev]) {
            handleBlock(prev, i);
            prev = i;
        }
    }
    handleBlock(prev, m_nAtoms);
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
    for(size_t idx = l; idx < r; ++idx) {
        for(size_t jdx = idx + 1; jdx < r; ++jdx) {
            size_t i = m_indicies[idx];
            size_t j = m_indicies[jdx];
            num_t dst = 0;
            for(size_t d = 0; d < UniverseDim; ++d) {
                dst += (m_coords[d][i] - m_coords[d][j]) * (m_coords[d][i] - m_coords[d][j]);
            }
            if(dst < (m_radiuses[i] + m_radiuses[j]) * (m_radiuses[i] + m_radiuses[j])) {
                m_collisionList.push_back({i, j});
            }
        }
    }
}


}
