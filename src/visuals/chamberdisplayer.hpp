#ifndef UNIVERSEDISPLAYER_HPP
#define UNIVERSEDISPLAYER_HPP

#include "chamber.hpp"
#include <QWidget>
#include <QVector>

class QSpinBox;
class QPushButton;

class ChamberDisplayer : public QWidget
{
    Q_OBJECT
public:
    explicit ChamberDisplayer(phys::Chamber::Metrics& metrics, QWidget *parent = nullptr);
    ~ChamberDisplayer() override;

    [[nodiscard]] const phys::Chamber::Metrics& getUniverseMetrics() const;

    void setScale(phys::LengthVal scale) { m_scale = scale; }

    enum class ColorPolicy {
        SingleColor,
        HeatColor,
    };

    ColorPolicy m_colorPolicy = ColorPolicy::HeatColor;

private:
    phys::Chamber::Metrics& m_chamberMetrics;
    phys::LengthVal m_scale;
    QTimer* m_timer;

    void resizeEvent(QResizeEvent* event) override;

    void rescale();

    QColor getColor(const phys::GasAtom&) const;

public slots:
    void paintEvent(QPaintEvent* event) override;
};

#endif // UNIVERSEDISPLAYER_HPP
