#include "mainwindow.hpp"
#include "ui_mainwindow.h"

#include "physicsthread.hpp"
#include "physconstants.hpp"
#include <QTimer>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_chamber({1e-8_m, 1e-8_m, 1e-9_m}),
    m_physThread(new PhysicsThread(m_chamber, this))
{

    m_cd = new ChamberDisplayer(m_chamberMetrics, this);
    m_cd->setGeometry(rect());
    m_cd->setScale(1e-8_m);
//    m_chamber.fillRandom(400, 1e-7_m / 1_sec, phys::num_t{4} * phys::consts::Dalton, 31e-12_m);
    m_chamber.fillRandomAxis(1000, 1e3_m / 1_sec, phys::num_t{4} * phys::consts::Dalton, 31e-12_m);

    m_timer = new QTimer(this);
    m_timer->setInterval(1000 / 60); // 60 fps
    m_timer->setSingleShot(false);
    connect(m_timer,SIGNAL(timeout()), this, SLOT(updateMetrics()));
    m_timer->start();

    m_chamber.setDT(1e-14_sec);
    m_physThread->setPeriod(0);
    
    ui->setupUi(this);

    connect(ui->startButton, SIGNAL(toggled(bool)), this, SLOT(toggleSimulation(bool)));
    connect(ui->timerBox, SIGNAL(valueChanged(int)), this, SLOT(setSimulationSpeed(int)));

    m_eDisplays[0] = ui->eDisplay1;
    m_eDisplays[1] = ui->eDisplay2;
    m_eDisplays[2] = ui->eDisplay3;

    m_pDisplays[0] = ui->pDisplay11;
    m_pDisplays[1] = ui->pDisplay12;
    m_pDisplays[2] = ui->pDisplay21;
    m_pDisplays[3] = ui->pDisplay22;
    m_pDisplays[4] = ui->pDisplay31;
    m_pDisplays[5] = ui->pDisplay32;

//    m_physThread->cont();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::resizeEvent(QResizeEvent *event)
{
    QMainWindow::resizeEvent(event);
    m_cd->setGeometry(rect());
}

void MainWindow::toggleSimulation(bool run)
{
    if(run)
        m_physThread->cont();
    else
        m_physThread->stop();
}

void MainWindow::setSimulationSpeed(int x)
{
    m_physThread->setPeriod(x);
}

void MainWindow::updateMetrics() {
    m_physThread->acquireMetrics(m_chamberMetrics);
    m_cd->update();

    QString str;
    QTextStream ss(&str);

    phys::Energy totalE{};

    for(size_t i = 0; i < phys::UniverseDim; i++) {
        ss << m_chamberMetrics.kineticEnergy[i];
        totalE += m_chamberMetrics.kineticEnergy[i];
        m_eDisplays[i]->setText(str);
        str.clear();
    }

    ss << totalE;
    ui->eDisplayTotal->setText(str);
    str.clear();

    ss << totalE * (phys::num_t{2./3.} / phys::num_t{m_chamberMetrics.atoms.size()}) / phys::consts::k;
    ui->tempDIsplay->setText(str);
    str.clear();

    for(size_t i = 0; i < 2 * phys::UniverseDim; i++) {
        ss << m_chamberMetrics.pressure[i];
        m_pDisplays[i]->setText(str);
        str.clear();
    }
//    ss << m_ud->getUniverseMetrics().energy;
//     ui->eDisplay->setText(str);
}


