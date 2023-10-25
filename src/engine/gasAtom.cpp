#include "gasAtom.hpp"
#include <iostream>

namespace phys {

void GasAtom::move(Time dt) {
    // std::cerr << m_pos << '\n';
    m_pos += m_v * dt;
}

}
