//
// Created by Laura Balascuta on 07/03/2021.
//

#include "OrangeTriangleSolver.h"

#include "../Utils/IO.h"
#include "../Utils/Math.h"
#include "../Entities/TriangleFactory.h"

#include <stdexcept>
#include <vector>

namespace Solver {

const int OrangeTriangleSolver::NO_INPUT_TRIANGLES{2};
const int OrangeTriangleSolver::NO_COORD_PER_POINT{3};
const std::string OrangeTriangleSolver::COORDINATES_FILENAME{"input.txt"};

OrangeTriangleSolver::OrangeTriangleSolver() : m_A1B1C1(nullptr),
                                               m_A2B2C2(nullptr),
                                               m_factory(new Entities::TriangleFactory()),
                                               m_solved(false) {

    // read the input triangles' coordinates from the text file
    std::vector<double> allCoordinates = Utils::IO::readCoordinatesFromFile(COORDINATES_FILENAME);

    // check if there are enough coordinates and populate the input vectors
    if (allCoordinates.size() == Entities::Triangle::POINTS_NO * NO_COORD_PER_POINT * NO_INPUT_TRIANGLES) {
        size_t i = 0;
        std::vector<double> a1b1c1;
        std::vector<double> a2b2c2;
        for (auto coord: allCoordinates) {
            ++i;
            if (i <= Entities::Triangle::POINTS_NO * 3) {
                a1b1c1.push_back(coord);
            } else {
                a2b2c2.push_back(coord);
            }
        }

        // populate the point's names for the input triangles
        std::vector<std::string> a1b1c1Name{"A1","B1","C1"};
        std::vector<std::string> a2b2c2Name{"A2","B2","C2"};

        // create the input green and blue triangles
        m_A1B1C1 = (std::dynamic_pointer_cast<Entities::Triangle>)(m_factory->createEntity(a1b1c1, a1b1c1Name));
        m_A2B2C2 = (std::dynamic_pointer_cast<Entities::Triangle>)(m_factory->createEntity(a2b2c2, a2b2c2Name));
    } else {
        std::cout << "Solver: There are not enough coordinates provided to create the input triangles." << std::endl;
    }
}

void OrangeTriangleSolver::solve() {
    // check if the triangles were created
    if (!m_A1B1C1 || !m_A2B2C2) {
        std::cout << "Solver: The input triangles could not be built." << std::endl;
        return;
    }

    m_solved = false;

    printInput();

    // the output coordinates for the orange triangles
    // (there are 2 values for each coordinate because the quadratic equation can have up to 2 roots)
    std::optional<double> xA3_1, xA3_2, yA3_1, yA3_2;
    std::optional<double> xB3_1, xB3_2, yB3_1, yB3_2;
    std::optional<double> xC3_1, xC3_2, yC3_1, yC3_2;

    try {
        // see the documentation for the algorithm and the result parameters and formulas
        double a, b, c, d, e, f, g, h, i, j;

        // calculate a and b
        ab(a, b);

        // calculate c, d, e
        cde(a, b, c, d, e);

        // calculate yA3_1 and yA3_2
        yA3(c, d, e, yA3_1, yA3_2);

        // calculate xA3_1 and xA3_2
        xA3(yA3_1, yA3_2, a, b, xA3_1, xA3_2);

        // calculate xB3_1 and xB3_2
        xB3(xA3_1, xA3_2, xB3_1, xB3_2);

        // calculate yB3_1 and yB3_2
        yB3(yA3_1, yA3_2, yB3_1, yB3_2);

        // calculate f and g
        fg(xA3_1, xB3_1, yA3_1, yB3_1, f, g);

        // calculate h, i and j
        hij(xA3_1, yA3_1, f, g, h, i, j);

        // calculate yC3_1 and yC3_2
        yC3(h, i, j, yC3_1, yC3_2);

        // calculate xC3_1 and xC3_2
        xC3(yC3_1, yC3_2, f, g, xC3_1, xC3_2);

        // mark it as solved if the algorithm is done
        m_solved = true;

    } catch (const std::exception &ex) {
        std::cout << "Solver: Exception occurred" << std::endl << ex.what() << std::endl;
    } catch (...) {
        std::cout << "Solver: Undefined exception occurred" << std::endl << std::endl;
    }

    // check for at least one output triangle
    if (!xA3_1.has_value() || !yA3_1.has_value()
        || !xB3_1.has_value() || !yB3_1.has_value()
        || !xC3_1.has_value() || !yC3_1.has_value()) {
        m_solved = false;
    }

    if (m_solved) {
        // combine the solutions and print out the orange triangle/s
        std::vector<std::string> a3b3c3Name{"A3","B3","C3"};
        std::vector<double> a3b3c3_1{xA3_1.value(), yA3_1.value(), m_A2B2C2->getABVector().getStartPoint().getZ(), //A
                                     xB3_1.value(), yB3_1.value(), m_A2B2C2->getABVector().getEndPoint().getZ(),   //B
                                     xC3_1.value(), yC3_1.value(), m_A2B2C2->getACVector().getEndPoint().getZ()};  //C

        addTriangleToOutputList(a3b3c3_1, a3b3c3Name);

        if (xC3_2.has_value() && yC3_2.has_value()) {
            std::vector<double>
                a3b3c3_2{xA3_1.value(), yA3_1.value(), m_A2B2C2->getABVector().getStartPoint().getZ(), //A
                         xB3_1.value(), yB3_1.value(), m_A2B2C2->getABVector().getEndPoint().getZ(),   //B
                         xC3_2.value(), yC3_2.value(), m_A2B2C2->getACVector().getEndPoint().getZ()};  //C

            addTriangleToOutputList(a3b3c3_2, a3b3c3Name);
        }

        if (xA3_2.has_value() && yA3_2.has_value() && xB3_2.has_value() && yB3_2.has_value()) {
            std::vector<double>
                a3b3c3_3{xA3_2.value(), yA3_2.value(), m_A2B2C2->getABVector().getStartPoint().getZ(), //A
                         xB3_2.value(), yB3_2.value(), m_A2B2C2->getABVector().getEndPoint().getZ(),   //B
                         xC3_1.value(), yC3_1.value(), m_A2B2C2->getACVector().getEndPoint().getZ()};  //C

            addTriangleToOutputList(a3b3c3_3, a3b3c3Name);

            if (xC3_2.has_value() && yC3_2.has_value()) {
                std::vector<double>
                    a3b3c3_4{xA3_2.value(), yA3_2.value(), m_A2B2C2->getABVector().getStartPoint().getZ(), //A
                             xB3_2.value(), yB3_2.value(), m_A2B2C2->getABVector().getEndPoint().getZ(),   //B
                             xC3_2.value(), yC3_2.value(), m_A2B2C2->getACVector().getEndPoint().getZ()};  //C

                addTriangleToOutputList(a3b3c3_4, a3b3c3Name);
            }
        }
    }

    printOutput();
}

void OrangeTriangleSolver::ab(double &a, double &b) {
    double denominator1 = m_A1B1C1->getABVector().getPositionVectorPoint().getX();
    double denominator2 = m_A1B1C1->getABVector().getPositionVectorPoint().getX();

    if (denominator1 == 0 || denominator2 == 0) {
        throw std::runtime_error("Attempted to divide by zero\n");
    }

    a = m_A1B1C1->getABVector().getPositionVectorPoint().getY() / denominator1;
    b = (m_A2B2C2->getABVector().getMiddlePoint().getX() * m_A1B1C1->getABVector().getPositionVectorPoint().getX() +
        m_A2B2C2->getABVector().getMiddlePoint().getY() * m_A1B1C1->getABVector().getPositionVectorPoint().getY() +
        m_A2B2C2->getABVector().getMiddlePoint().getZ() * m_A1B1C1->getABVector().getPositionVectorPoint().getZ() -
        m_A2B2C2->getABVector().getStartPoint().getZ() * m_A1B1C1->getABVector().getPositionVectorPoint().getZ()) /
        denominator2;
}

void OrangeTriangleSolver::cde(const double &a, const double &b, double &c, double &d, double &e) {
    c = pow(a, 2) + 1;
    d = 2 * m_A2B2C2->getABVector().getMiddlePoint().getX() * a
        - 2 * a * b
        - 2 * m_A2B2C2->getABVector().getMiddlePoint().getY();
    e = -2 * m_A2B2C2->getABVector().getMiddlePoint().getX() * b + pow(b, 2)
        - pow(m_A1B1C1->getABVector().getMagnitude(), 2) / 4
        + pow(m_A2B2C2->getABVector().getEndPoint().getZ() - m_A2B2C2->getABVector().getStartPoint().getZ(), 2) / 4
        + pow(m_A2B2C2->getABVector().getMiddlePoint().getX(), 2)
        + pow(m_A2B2C2->getABVector().getMiddlePoint().getY(), 2);
}

void OrangeTriangleSolver::yA3(const double &c,
                               const double &d,
                               const double &e,
                               std::optional<double> &yA3_1,
                               std::optional<double> &yA3_2) {
    Utils::Math::quadraticEquation(c, d, e, yA3_1, yA3_2);
}

void OrangeTriangleSolver::xA3(const std::optional<double> &yA3_1, const std::optional<double> &yA3_2,
                               const double &a, const double &b,
                               std::optional<double> &xA3_1, std::optional<double> &xA3_2) {
    if (yA3_1.has_value()) {
        xA3_1 = -yA3_1.value() * a + b;
    }

    if (yA3_2.has_value()) {
        xA3_2 = -yA3_2.value() * a + b;
    }
}

void OrangeTriangleSolver::xB3(const std::optional<double> &xA3_1, const std::optional<double> &xA3_2,
                               std::optional<double> &xB3_1, std::optional<double> &xB3_2) {
    if (xA3_1.has_value()) {
        xB3_1 = 2 * m_A2B2C2->getABVector().getMiddlePoint().getX() - xA3_1.value();
    }

    if (xA3_2.has_value()) {
        xB3_2 = 2 * m_A2B2C2->getABVector().getMiddlePoint().getX() - xA3_2.value();
    }
}

void OrangeTriangleSolver::yB3(const std::optional<double> &yA3_1, const std::optional<double> &yA3_2,
                               std::optional<double> &yB3_1, std::optional<double> &yB3_2) {
    if (yA3_1.has_value()) {
        yB3_1 = 2 * m_A2B2C2->getABVector().getMiddlePoint().getY() - yA3_1.value();
    }

    if (yA3_2.has_value()) {
        yB3_2 = 2 * m_A2B2C2->getABVector().getMiddlePoint().getY() - yA3_2.value();
    }
}

void OrangeTriangleSolver::fg(const std::optional<double> &xA3, const std::optional<double> &xB3,
                              const std::optional<double> &yA3, const std::optional<double> &yB3,
                              double &f, double &g) {
    // there is no need to use both _1 and _2 roots because the value will be the same
    if (xA3.has_value() && xB3.has_value() && yA3.has_value() && yB3.has_value()) {
        double denominator1 = -xB3.value() + xA3.value();
        double denominator2 = -2 * xB3.value() + 2 * xA3.value();

        if (denominator1 == 0 || denominator2 == 0) {
            throw std::runtime_error("Attempted to divide by zero\n");
        }

        f = (yB3.value() - yA3.value()) / denominator1;
        g = pow(m_A1B1C1->getBCVector().getMagnitude(), 2)
            - pow(m_A2B2C2->getABVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getEndPoint().getZ(), 2)
            - pow(xB3.value(), 2) - pow(yB3.value(), 2)
            - pow(m_A1B1C1->getACVector().getMagnitude(), 2)
            + pow(m_A2B2C2->getACVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getStartPoint().getZ(), 2)
            + pow(xA3.value(), 2)
            + pow(yA3.value(), 2);
        g = g / denominator2;
    }
}

void OrangeTriangleSolver::hij(const std::optional<double> &xA3, const std::optional<double> &yA3,
                               const double &f, const double &g,
                               double &h, double &i, double &j) {
    if (xA3.has_value() && yA3.has_value()) {
        h = pow(f, 2) + 1;
        i = 2 * (f * g - xA3.value() * f - yA3.value());
        j = pow(g, 2) - 2 * xA3.value() * g
            - pow(m_A1B1C1->getACVector().getMagnitude(), 2)
            + pow(m_A2B2C2->getACVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getStartPoint().getZ(), 2)
            + pow(xA3.value(), 2) + pow(yA3.value(), 2);
    }
}

void OrangeTriangleSolver::yC3(const double &h, const double &i, const double &j,
                               std::optional<double> &yC3_1,
                               std::optional<double> &yC3_2) {
    Utils::Math::quadraticEquation(h, i, j, yC3_1, yC3_2);
}

void OrangeTriangleSolver::xC3(const std::optional<double> &yC3_1, const std::optional<double> &yC3_2,
                               const double &f, const double &g,
                               std::optional<double> &xC3_1, std::optional<double> &xC3_2) {
    if (yC3_1.has_value()) {
        xC3_1 = yC3_1.value() * f + g;
    }

    if (yC3_2.has_value()) {
        xC3_2 = yC3_2.value() * f + g;
    }
}

void OrangeTriangleSolver::addTriangleToOutputList(const std::vector<double> &coord, const std::vector<std::string> &name) {
    std::shared_ptr<Entities::Triangle> outputTriangle = (std::dynamic_pointer_cast<Entities::Triangle>)
        (m_factory->createEntity(coord, name));

    if (checkAssumptionsSolved(outputTriangle)) {
        m_outputTriangles.push_back(outputTriangle);
    }
}

bool OrangeTriangleSolver::checkAssumptionsSolved(const std::shared_ptr<Entities::Triangle> &triangle) const {
    return  checkSamePlane(triangle) &&
            checkMagnitudeEquality(triangle) &&
        checkPerpendicularVectors(triangle) &&
            checkMiddle(triangle);
}

bool OrangeTriangleSolver::checkSamePlane(const std::shared_ptr<Entities::Triangle> &triangle) const {
    return triangle->getABVector().getStartPoint().getZ() == m_A2B2C2->getABVector().getStartPoint().getZ() &&
           triangle->getABVector().getEndPoint().getZ() == m_A2B2C2->getABVector().getEndPoint().getZ() &&
           triangle->getBCVector().getEndPoint().getZ() == m_A2B2C2->getABVector().getEndPoint().getZ();
}

bool OrangeTriangleSolver::checkMagnitudeEquality(const std::shared_ptr<Entities::Triangle> &triangle) const {
    return triangle->getABVector().getMagnitude() == m_A1B1C1->getABVector().getMagnitude() &&
           triangle->getBCVector().getMagnitude() == m_A1B1C1->getBCVector().getMagnitude() &&
           triangle->getACVector().getMagnitude() == m_A1B1C1->getACVector().getMagnitude();
}

bool OrangeTriangleSolver::checkPerpendicularVectors(const std::shared_ptr<Entities::Triangle> &triangle) const {
    return triangle->getABVector().computeDotProduct(m_A1B1C1->getABVector()) == 0;
}

bool OrangeTriangleSolver::checkMiddle(const std::shared_ptr<Entities::Triangle> &triangle) const {
    return triangle->getABVector().getMiddlePoint().getX() == m_A2B2C2->getABVector().getMiddlePoint().getX() &&
           triangle->getABVector().getMiddlePoint().getY() == m_A2B2C2->getABVector().getMiddlePoint().getY() &&
           triangle->getABVector().getMiddlePoint().getZ() == m_A2B2C2->getABVector().getMiddlePoint().getZ();
}

void OrangeTriangleSolver::printInput() const {
    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << std::endl;
    std::cout << "INPUT TRIANGLES" << std::endl;
    std::cout << std::endl;
    std::cout << "1. green triangle:" << std::endl;
    m_A1B1C1->printInfo();
    std::cout << std::endl;
    std::cout << "2. blue triangle:" << std::endl;
    m_A2B2C2->printInfo();
    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << std::endl;
}

void OrangeTriangleSolver::printOutput() const {
    std::cout << "OUTPUT" << std::endl;

    if (m_solved) {
        std::cout << std::endl;
        std::cout << "orange triangle/s:" << std::endl;

        std::size_t i=0;
        for (const auto& triangle: m_outputTriangles) {
            std::cout << ++i << ". ";
            triangle->printInfo();
        }

    } else {
        std::cout << std::endl;
        std::cout << "Solver: The problem has no solutions!" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "---------------" << std::endl;
    std::cout << std::endl;
}

}