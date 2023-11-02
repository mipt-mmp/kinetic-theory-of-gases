#ifndef ENGINE_BALLSCOLLECTION_HPP
#define ENGINE_BALLSCOLLECTION_HPP
#include "gasAtom.hpp"
#include "units.hpp"

#include <functional>
#include <QMutex>

namespace phys {

class BallsCollection;

namespace detail {
    class GasAtomProxy {
        BallsCollection& m_balls;
        size_t m_index;
        GasAtom m_atom;

        GasAtomProxy(BallsCollection& balls, size_t index);

        void updateBalls();

        friend BallsCollection;
    public:
        const Mass& getMass() const {
            return m_atom.getMass();
        }

        Length getRadius() const {
            return m_atom.getRadius();
        }

        void collide(Time t) {
            m_atom.collide(t);
            updateBalls();
        }

        Length getFreeFlight(Time t) {
            return m_atom.getFreeFlight(t);
        }

        Energy getAverageEnergy(Time t) {
            return m_atom.getAverageEnergy(t);
        }

        const Velocity& getVelocity() const {
            return m_atom.getVelocity();
        }

        void setVelocity(const Velocity& newVelocity) {
            m_atom.setVelocity(newVelocity);
            updateBalls();
        }

        const Position& getPos() const {
            return m_atom.getPos();
        }

        void setPos(const Position& newPos) {
            m_atom.setPos(newPos);
            updateBalls();
        }

        void move(Time dt) {
           m_atom.move(dt);
           updateBalls();
        }

        ImpulseMoment getImpulseMoment() const {
            return m_atom.getImpulseMoment();
        }

        Impulse getImpulse() const {
            return m_atom.getImpulse();
        }

        Energy getKinetic() const {
           return m_atom.getKinetic();
        }

        Vector<Energy> getKineticDistributed() const {
            return m_atom.getKineticDistributed();
        }

    };
}

class BallsCollection {
    std::array<std::vector<num_t>, UniverseDim> m_coords;
    std::array<std::vector<num_t>, UniverseDim> m_velocities;
    std::vector<num_t> m_masses;
    std::vector<num_t> m_radiuses;
    size_t m_nAtoms = 0;


    friend detail::GasAtomProxy;

    Length m_mScale;
    Time   m_tScale;

    std::array<num_t, UniverseDim> m_walls;
public:
    static const std::size_t MeasurementSize = 64;
private:
    std::array<std::array<num_t, 2 * UniverseDim>, MeasurementSize> m_wallImpulse;
    std::size_t m_stepIdx = 0;

    std::vector<uint32_t> m_hashes;
    std::vector<uint32_t> m_indicies;
    std::vector<uint32_t> m_radixBuffer;
    std::vector<uint32_t> m_radixIndiciesBuffer;

    std::vector<uint32_t> m_cellCounter;

    num_t m_cellSize;
    std::array<uint32_t, UniverseDim> m_shifts;

    std::vector<std::pair<size_t, size_t>> m_collisionList;
    QMutex m_listMutex;

    bool m_enableHole = false;

public:
    BallsCollection(Length meterScale, Time timeScale) : m_mScale(meterScale), m_tScale(timeScale) {}

    ImpulseVal getWallImpulse(size_t i) const {
        num_t val;
        for(size_t t = 0; t < MeasurementSize;t++) {
            val += m_wallImpulse[t][i];
        }
        return std::abs(val) * m_mScale / m_tScale * Mass{1};
    }

    template<typename F>
    void addAtoms(size_t N, F generator) {
        for(size_t i = 0; i < N; ++i) {
            GasAtom atom = generator();
            for(size_t j = 0; j < UniverseDim; ++j) {
                m_coords    [j].push_back(*(atom.getPos()     [j] / m_mScale));
                m_velocities[j].push_back(*(atom.getVelocity()[j] / m_mScale * m_tScale));
            }
            m_masses  .push_back(*atom.getMass());
            m_radiuses.push_back(*(atom.getRadius() / m_mScale));
            m_nAtoms++;
        }

        m_hashes             .resize(m_nAtoms);
        m_indicies           .resize(m_nAtoms);
        m_radixBuffer        .resize(m_nAtoms);
        m_radixIndiciesBuffer.resize(m_nAtoms);
    }

    void setWalls(Position pos) {
        for(size_t i = 0; i < UniverseDim; ++i) {
            m_walls[i] = *(pos[i] / m_mScale);
        }
    }

    void setCellSize(Length l);

    detail::GasAtomProxy operator[](size_t i) {return detail::GasAtomProxy(*this, i);}

    void deleteAtom(size_t i);

    GasAtom getAtom(size_t i) const;

    void move(Time dt);

    size_t size() const {return m_nAtoms;}

    void push_back(const GasAtom& atom) {
        addAtoms(1, [&](){return atom;});
    }

    void handleWallCollisions();
    
    void handleCollisions();

    const std::vector<std::pair<size_t, size_t>>& getCollisions() { return m_collisionList; }
    void setEnableHole(bool newEnableHole);

private:
    void radixSort();
    
    void handleSub(size_t i, size_t j);

    void handleBlock(size_t i, size_t j);
};

}

#endif /* ENGINE_BALLSCOLLECTION_HPP */
