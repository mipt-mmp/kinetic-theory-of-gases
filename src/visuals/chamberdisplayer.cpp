#include "chamberdisplayer.hpp"
#include <QTimer>
#include <QSpinBox>
#include <QPushButton>
#include <QDateTimeEdit>
#include <QPainter>
#include <QDebug>
#include <QBrush>
#include <QPen>
#include <iostream>
ChamberDisplayer::ChamberDisplayer(phys::Chamber::Metrics& metrics, QWidget *parent) : QWidget(parent), m_chamberMetrics(metrics) {
    rescale();

    QPalette pal = QPalette();
    pal.setColor(QPalette::Window, Qt::white);
    setAutoFillBackground(true);
    setPalette(pal);
}

ChamberDisplayer::~ChamberDisplayer() {}

void ChamberDisplayer::resizeEvent(QResizeEvent* event)
{
    QWidget::resizeEvent(event);
    rescale();
}

void ChamberDisplayer::rescale()
{}

void ChamberDisplayer::paintEvent(QPaintEvent* /*event*/)
{
    auto& atoms = m_chamberMetrics.atoms;
    if constexpr (phys::UniverseDim >= 3) {
        // std::sort(atoms.begin(), atoms.end(), [](const phys::GasAtom& lhs, const phys::GasAtom& rhs) -> bool{
        //     return lhs.getPos().Z() < rhs.getPos().Z();
        // });
    } 
    QPainter painter(this);
    QPen pen;
    pen.setWidth(3);
    painter.setPen(pen);
    phys::num_t pixscale{std::min(rect().width(), rect().height())};
    
    painter.drawRoundedRect(0, 0, 
        static_cast<int>(pixscale * *(m_chamberMetrics.chamberCorner.X() / m_scale)),
        static_cast<int>(pixscale * *(m_chamberMetrics.chamberCorner.Y() / m_scale)),
        3, 3);

    painter.setPen(pen);
    QBrush brush(Qt::SolidPattern);
    painter.setBrush(brush);
    for(auto& atom : atoms) {
        QColor color = getColor(atom);
        brush.setColor(color);
        pen.setColor  (color);
        painter.setPen(pen);
        painter.setBrush(brush);
        
        int radius = static_cast<int>(pixscale * *(atom.getRadius()  / m_scale));
        painter.drawEllipse(
            static_cast<int>(pixscale * *((atom.getPos().X()) / m_scale)) - radius,
            static_cast<int>(pixscale * *((atom.getPos().Y()) / m_scale)) - radius,
            2 * radius,
            2 * radius);
    }
}

QColor ChamberDisplayer::getColor(const phys::GasAtom& atom) const {
    switch (m_colorPolicy) {

    case ColorPolicy::SingleColor:
        return Qt::gray;
    case ColorPolicy::HeatColor:
    {
        int hue = std::min(200, static_cast<int>(*atom.getKinetic() * 300'000));
        return QColor::fromHsv(200 - hue, 250, 250);
    }
    }
    return Qt::magenta;
}
