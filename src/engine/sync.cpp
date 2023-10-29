#include "sync.hpp"

namespace QSync {

QCyclicBarrier::QCyclicBarrier(size_t count)
    : m_maxCount(count)
    , m_curCount(0)
    , m_generation(0) {}

void QCyclicBarrier::wait() {
    m_mutex.lock();
    size_t thisGeneration = m_generation;

    if (++m_curCount == m_maxCount) {
        m_generation++;
        m_curCount = 0;
        m_condvar.wakeAll();
    } else {
        while (thisGeneration == m_generation) {
            m_condvar.wait(&m_mutex);
        }
    }

    m_mutex.unlock();
}

QBarrier::QBarrier(int32_t count)
    : m_impl(new QCyclicBarrier(count)) {}

void QBarrier::wait() {
    m_impl->wait();
}

} // namespace QSync
