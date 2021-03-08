//
// Created by Laura Balascuta on 07/03/2021.
//

#include "Triangle.h"

namespace Entities {

const int Triangle::POINTS_NO = 3;

Triangle::Triangle(const Vector& AB, const Vector& BC, const Vector& AC) :
    m_AB(AB),
    m_BC(BC),
    m_AC(AC),
    m_name(AB.getStartPoint().getName() + AB.getEndPoint().getName() + AC.getEndPoint().getName()) {
}

Triangle::~Triangle() {
}

void Triangle::printInfo() const {
    std::cout << *this << std::endl;
}

std::string Triangle::getName() const {
    return m_name;
}

Vector Triangle::getABVector() const {
    return m_AB;
}

Vector Triangle::getBCVector() const {
    return m_BC;
}

Vector Triangle::getACVector() const {
    return m_AC;
}

std::ostream &operator<<(std::ostream &os, const Triangle &triangle) {
    os << triangle.getName() <<"[" << triangle.getABVector().getStartPoint() << ", " << triangle.getBCVector().getStartPoint() << ", " << triangle.getACVector().getEndPoint() << "]";
    return os;
}

}