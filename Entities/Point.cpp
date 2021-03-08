//
// Created by Laura Balascuta on 07/03/2021.
//

#include "Point.h"

namespace Entities {

Point::Point(const double &x, const double &y, const double &z, const std::string &name) : m_x(x), m_y(y), m_z(z),
                                                                                           m_name(name) {
}

void Point::printInfo() const {
    std::cout << *this << std::endl;
}

std::string Point::getName() const {
    return m_name;
}

double Point::getX() const {
    return m_x;
}

double Point::getY() const {
    return m_y;
}

double Point::getZ() const {
    return m_z;
}

std::ostream &operator<<(std::ostream &os, const Point &point) {
    os << point.getName() << "[" << point.getX() << ", " << point.getY() << ", " << point.getZ() << "]";
    return os;
}

}