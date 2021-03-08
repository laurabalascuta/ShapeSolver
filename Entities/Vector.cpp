//
// Created by Laura Balascuta on 07/03/2021.
//

#include "Vector.h"

#include <math.h>

namespace Entities {

const std::string Vector::MIDDLE_POINT_NAME = "M";

Vector::Vector(const Point &start, const Point &end) :
    m_start(start),
    m_end(end),
    m_positionVector(end.getX() - start.getX(),
                     end.getY() - start.getY(),
                     end.getZ() - start.getZ(),
                     start.getName() + end.getName()),
    m_middle((start.getX() + end.getX()) / 2,
             (start.getY() + end.getY()) / 2,
             (start.getZ() + end.getZ()) / 2,
             MIDDLE_POINT_NAME),
    m_name(start.getName() + end.getName()) {
    m_magnitude = sqrt(pow((m_end.getX() - m_start.getX()), 2) +
                       pow((m_end.getY() - m_start.getY()), 2) +
                       pow((m_end.getZ() - m_start.getZ()), 2));
}

void Vector::printInfo() const {
    std::cout << *this << std::endl;
}

std::string Vector::getName() const {
    return m_name;
}

Point Vector::getStartPoint() const {
    return m_start;
}

Point Vector::getEndPoint() const {
    return m_end;
}

Point Vector::getPositionVectorPoint() const {
    return m_positionVector;
}

Point Vector::getMiddlePoint() const {
    return m_middle;
}

double Vector::getMagnitude() const {
    return m_magnitude;
}

std::ostream &operator<<(std::ostream &os, const Vector &vector) {
    os << vector.getName() << "[" << vector.getStartPoint() << ", " << vector.getEndPoint() << "]";
    return os;
}

}