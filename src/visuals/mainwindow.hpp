#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <array>
#include "chamberdisplayer.hpp"

namespace Ui {
class MainWindow;
}
class PhysicsThread;
class QTimer;
class QLineEdit;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

    void resizeEvent(QResizeEvent* event) override;

private:
    ChamberDisplayer* m_cd;
    Ui::MainWindow *ui;
    
    QTimer* m_timer;

    std::array<QLineEdit*, 3> m_eDisplays;
    std::array<QLineEdit*, 6> m_pDisplays;
    
    phys::Chamber m_chamber;
    phys::Chamber::Metrics m_chamberMetrics;

    PhysicsThread* m_physThread;

private slots:
    void toggleSimulation(bool);
    void setSimulationSpeed(int);

    void updateMetrics();
};

#endif // MAINWINDOW_H
