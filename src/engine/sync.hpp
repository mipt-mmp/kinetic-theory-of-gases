#ifndef ENGINE_SYNC_HPP
#define ENGINE_SYNC_HPP

#include <QMutex>
#include <QWaitCondition>
#include <QSharedPointer>

namespace QSync {

class QCyclicBarrier {
public:
    explicit QCyclicBarrier(size_t count);

    void wait();

private:
    Q_DISABLE_COPY(QCyclicBarrier)

    size_t m_maxCount;
    size_t m_curCount;
    size_t m_generation;

    QMutex m_mutex;
    QWaitCondition m_condvar;
};

class QBarrier {
public:
    explicit QBarrier(int32_t count = 0);

    void wait();

private:
    QSharedPointer<QCyclicBarrier> m_impl;
};

} // namespace QSync

#endif // ENGINE_SYNC_HPP
