#include "mainwindow.hpp"
#include "ui_mainwindow.h"

#include "physicsthread.hpp"

#include <QTimer>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_chamber({1_m, 1_m}),
    m_physThread(new PhysicsThread(m_chamber, this))
{
    ui->setupUi(this);

    m_cd = new ChamberDisplayer(m_chamberMetrics, this);
    m_cd->setGeometry(rect());
    m_cd->setScale(1_m);
    m_chamber.fillRandom(300, 0.1_m / 1_sec, 1_kg, 0.002_m);

    m_timer = new QTimer(this);
    m_timer->setInterval(1000 / 60); // 60 fps
    m_timer->setSingleShot(false);
    connect(m_timer,SIGNAL(timeout()), this, SLOT(updateMetrics()));
    m_timer->start();

    m_chamber.setDT(0.001_sec);
    m_physThread->setPeriod(0);
    
    m_physThread->cont();
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

void MainWindow::updateMetrics() {
    m_physThread->acquireMetrics(m_chamberMetrics);
    m_cd->update();
}


