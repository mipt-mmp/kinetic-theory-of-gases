#include "physicsthread.hpp"
#include <QMutexLocker>
#include <iostream>

PhysicsThread::PhysicsThread(phys::Chamber& chamber, QObject* parent)
    : QThread(parent)
    , m_chamber(chamber) {
    start(LowPriority);
}

PhysicsThread::~PhysicsThread() {}

void PhysicsThread::run() {
    forever {
        {
            QMutexLocker lock(&m_mutex);

            if (m_stopped) {
                m_allow_run.wait(&m_mutex);
            }

            m_chamber.step();
        }

        if (m_period != -1) {
            long time = static_cast<long>(std::pow(1.2, m_period));
            time = std::min(time, 5'000'000l);
            if (time < 0)
                time = 5'000'000l;
            usleep(time);
        } else {
            for (size_t i = 0; i < 100'000; ++i) {
                __asm__ volatile("nop");
            }
        }
    }
}

[[nodiscard]] bool PhysicsThread::getStopped() {
    QMutexLocker lock(&m_mutex);
    return m_stopped;
}

int PhysicsThread::getPeriod() const {
    return m_period;
}

void PhysicsThread::setPeriod(int newPeriod) {
    QMutexLocker lock(&m_mutex);
    m_period = newPeriod;
}
