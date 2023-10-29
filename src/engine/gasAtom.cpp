#include "gasAtom.hpp"
#include <iostream>

namespace phys {

uint64_t GasAtom::getBlockId() const {
    return m_blockId;
}

void GasAtom::setBlockId(uint64_t blockId) {
    m_blockId = blockId;
}

bool GasAtom::operator<(const GasAtom& other) {
    return m_blockId < other.m_blockId;
}

} // namespace phys
