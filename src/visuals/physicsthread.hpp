#ifndef PHYSICSTHREAD_HPP
#define PHYSICSTHREAD_HPP

#include <QThread>
#include <QMutex>
#include <QWaitCondition>

#include "chamber.hpp"

class PhysicsThread : public QThread
{
    Q_OBJECT
public:
    PhysicsThread(phys::Chamber& universe, QObject* parent = nullptr);
    ~PhysicsThread() override;

public slots:
    void stop() {
        m_mutex.lock();
        m_stopped = true;
        m_mutex.unlock();
        emit toggled(false);
    }

    void cont() {
        m_mutex.lock();
        if(m_stopped) {
            m_stopped = false;
            m_allow_run.wakeAll();
            emit toggled(true);
        }
        m_mutex.unlock();
    }

    void toggle() {
        m_mutex.lock();
        m_stopped ^= 1;
        if(m_stopped) {
            m_mutex.unlock();
            emit toggled(false);
            return;
        }

        m_allow_run.wakeAll();
        emit toggled(true);
        m_mutex.unlock();
    }

public:
    [[nodiscard]] int getPeriod() const;
    [[nodiscard]] bool getStopped();

    void acquireMetrics(phys::Chamber::Metrics& metrics) {
        QMutexLocker lock(&m_mutex);
        m_chamber.getMetrics(metrics);
    }

signals:
    void toggled(bool);

public slots:
    void setPeriod(int newPeriod);

protected:
    void run() override;

private:
    phys::Chamber& m_chamber;
    QMutex m_mutex;
    QWaitCondition m_allow_run;
    bool m_stopped = true;
    int m_period = 32;
};

#endif // PHYSICSTHREAD_HPP
