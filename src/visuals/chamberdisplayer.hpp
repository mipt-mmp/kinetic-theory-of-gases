#ifndef UNIVERSEDISPLAYER_HPP
#define UNIVERSEDISPLAYER_HPP

#include "chamber.hpp"
#include <QVector>
#include <QWidget>

class QSpinBox;
class QPushButton;

class ChamberDisplayer : public QWidget {
    Q_OBJECT
public:
    explicit ChamberDisplayer(phys::Chamber::Metrics& metrics, QWidget* parent = nullptr);
    ~ChamberDisplayer() override;

    [[nodiscard]] const phys::Chamber::Metrics& getUniverseMetrics() const;

    void setScale(phys::LengthVal scale) {
        m_scale = scale;
    }

    enum class ColorPolicy {
        SingleColor,
        HeatColor,
    };

    ColorPolicy m_colorPolicy = ColorPolicy::HeatColor;


private:
    std::array<std::pair<QPoint, QColor>, 1024> m_record;
    std::size_t m_recordIdx = 0;
    std::size_t m_followIdx = 0;

    phys::Chamber::Metrics& m_chamberMetrics;
    phys::LengthVal m_scale;
    QTimer* m_timer;

    bool m_follow = false;

    void resizeEvent(QResizeEvent* event) override;

    void rescale();

    QColor getColor(const phys::GasAtom&) const;

public slots:
    void paintEvent(QPaintEvent* event) override;

    void setFollowIdx(int newFollowIdx);

    void setFollow(bool newFollow);
};

#endif // UNIVERSEDISPLAYER_HPP
