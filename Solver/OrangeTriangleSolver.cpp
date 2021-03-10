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
    std::vector<std::vector<double>> outputTrianglesCoordinates;
    try {
        double zA3 = m_A2B2C2->getABVector().getStartPoint().getZ();
        double zB3 = m_A2B2C2->getABVector().getEndPoint().getZ();
        double zC3 = m_A2B2C2->getACVector().getEndPoint().getZ();

        // see the documentation for the algorithm, the result parameters and formulas
        double a, b, c, d, e, f, g, h, i, j, k, l;

        if (m_A1B1C1->getABVector().getPositionVectorPoint().getX() == 0) {
            double yA3;
            compute_yA3(yA3);

            double m, n;
            compute_mn(yA3, m, n);

            std::vector<double> xA3_roots;
            compute_quadratic_xA3(m, n, xA3_roots);

            for (auto xA3: xA3_roots) {
                double xB3, yB3;
                compute_xB3(xA3, xB3);
                compute_yB3(yA3, yB3);

                if (xA3 != xB3) {
                    compute_hi(xA3, xB3, yA3, yB3, h, i);
                    // calculate j, k and l
                    compute_jkl(xA3, yA3, h, i, j, k, l);

                    std::vector<double> yC3_roots;
                    compute_quadratic_yC3(j, k, l, yC3_roots);
                    for (auto yC3: yC3_roots) {
                        double xC3;
                        compute_xC3(yC3, h, i, xC3);
                        std::vector<double> coord{xA3, yA3, zA3,
                                                  xB3, yB3, zB3,
                                                  xC3, yC3, zC3};
                        outputTrianglesCoordinates.push_back(coord);
                    }
                } else {
                    double yC3;
                    compute_yC3(xA3, xB3, yA3, yB3, yC3);

                    compute_fg(xA3, yA3, yC3, f, g);

                    std::vector<double> xC3_roots;
                    compute_quadratic_xC3(1, f, g, xC3_roots);
                    for (auto xC3: xC3_roots) {
                        std::vector<double> coord{xA3, yA3, zA3,
                                                  xB3, yB3, zB3,
                                                  xC3, yC3, zC3};
                        outputTrianglesCoordinates.push_back(coord);
                    }
                }
            }
        } else {

            compute_ab(a, b);

            compute_cde(a, b, c, d, e);

            std::vector<double> yA3_roots;
            compute_quadratic_yA3(c, d, e, yA3_roots);

            for (auto yA3: yA3_roots) {
                double xA3;
                compute_xA3(yA3, a, b, xA3);

                double xB3, yB3;
                compute_xB3(xA3, xB3);
                compute_yB3(yA3, yB3);

                if (xA3 != xB3) {
                    compute_hi(xA3, xB3, yA3, yB3, h, i);
                    compute_jkl(xA3, yA3, h, i, j, k, l);

                    std::vector<double> yC3_roots;
                    compute_quadratic_yC3(j, k, l, yC3_roots);
                    for (auto yC3: yC3_roots) {
                        double xC3;
                        compute_xC3(yC3, h, i, xC3);
                        std::vector<double> coord{xA3, yA3, zA3,
                                                  xB3, yB3, zB3,
                                                  xC3, yC3, zC3};
                        outputTrianglesCoordinates.push_back(coord);
                    }
                } else {
                    double yC3;
                    compute_yC3(xA3, xB3, yA3, yB3, yC3);

                    compute_fg(xA3, yA3, yC3, f, g);

                    std::vector<double> xC3_roots;
                    compute_quadratic_xC3(1, f, g, xC3_roots);
                    for (auto xC3: xC3_roots) {
                        std::vector<double> coord{xA3, yA3, zA3,
                                                  xB3, yB3, zB3,
                                                  xC3, yC3, zC3};
                        outputTrianglesCoordinates.push_back(coord);
                    }
                }
            }
        }

        // mark it as solved if the algorithm is done
        m_solved = true;

    } catch (const std::exception &ex) {
        std::cout << "Solver: Exception occurred" << std::endl << ex.what() << std::endl;
    } catch (...) {
        std::cout << "Solver: Undefined exception occurred" << std::endl << std::endl;
    }


    if (m_solved) {
        // combine the solutions and print out the orange triangle/s
        std::vector<std::string> a3b3c3Name{"A3","B3","C3"};
        for (const auto &coord: outputTrianglesCoordinates) {
            addTriangleToOutputList(coord, a3b3c3Name);
        }
    }

    printOutput();
}

void OrangeTriangleSolver::compute_ab(double &a, double &b) {
    double denominator = m_A1B1C1->getABVector().getPositionVectorPoint().getX();

    a = m_A1B1C1->getABVector().getPositionVectorPoint().getY() / denominator;
    b = (m_A2B2C2->getABVector().getMiddlePoint().getX() * m_A1B1C1->getABVector().getPositionVectorPoint().getX() +
        m_A2B2C2->getABVector().getMiddlePoint().getY() * m_A1B1C1->getABVector().getPositionVectorPoint().getY() +
        m_A2B2C2->getABVector().getMiddlePoint().getZ() * m_A1B1C1->getABVector().getPositionVectorPoint().getZ() -
        m_A2B2C2->getABVector().getStartPoint().getZ() * m_A1B1C1->getABVector().getPositionVectorPoint().getZ()) /
        denominator;
}

void OrangeTriangleSolver::compute_cde(const double &a, const double &b, double &c, double &d, double &e) {
    c = pow(a, 2) + 1;
    d = 2 * m_A2B2C2->getABVector().getMiddlePoint().getX() * a
        - 2 * a * b
        - 2 * m_A2B2C2->getABVector().getMiddlePoint().getY();
    e = - 2 * m_A2B2C2->getABVector().getMiddlePoint().getX() * b + pow(b, 2)
        - pow(m_A1B1C1->getABVector().getMagnitude(), 2) / 4
        + pow(m_A2B2C2->getABVector().getEndPoint().getZ() - m_A2B2C2->getABVector().getStartPoint().getZ(), 2) / 4
        + pow(m_A2B2C2->getABVector().getMiddlePoint().getX(), 2)
        + pow(m_A2B2C2->getABVector().getMiddlePoint().getY(), 2);
}

void OrangeTriangleSolver::compute_yA3(double &yA3) {
    yA3 = m_A2B2C2->getABVector().getMiddlePoint().getX() * m_A1B1C1->getABVector().getPositionVectorPoint().getX() +
          m_A2B2C2->getABVector().getMiddlePoint().getY() * m_A1B1C1->getABVector().getPositionVectorPoint().getY() +
          m_A2B2C2->getABVector().getMiddlePoint().getZ() * m_A1B1C1->getABVector().getPositionVectorPoint().getZ() -
          m_A2B2C2->getABVector().getStartPoint().getZ() * m_A1B1C1->getABVector().getPositionVectorPoint().getZ();
    yA3 = yA3 / m_A1B1C1->getABVector().getPositionVectorPoint().getY();
}

void OrangeTriangleSolver::compute_mn(const double &yA3, double &m, double &n) {
    m = - 2 *  m_A2B2C2->getABVector().getMiddlePoint().getX();
    n = - 2 * m_A2B2C2->getABVector().getMiddlePoint().getY() * yA3 + pow(yA3, 2)
        - pow(m_A1B1C1->getABVector().getMagnitude(), 2) / 4
        + pow(m_A2B2C2->getABVector().getEndPoint().getZ() - m_A2B2C2->getABVector().getStartPoint().getZ(), 2) / 4
        + pow(m_A2B2C2->getABVector().getMiddlePoint().getX(), 2)
        + pow(m_A2B2C2->getABVector().getMiddlePoint().getY(), 2);
}

void OrangeTriangleSolver::compute_quadratic_xA3(const double &m, const double &n, std::vector<double> &xA3) {
    Utils::Math::quadraticEquation(1, m, n, xA3);
}

void OrangeTriangleSolver::compute_quadratic_yA3(const double &c, const double &d, const double &e,
                                                 std::vector<double> &yA3) {
    Utils::Math::quadraticEquation(c, d, e, yA3);
}

void OrangeTriangleSolver::compute_xA3(const double &yA3, const double &a, const double &b, double &xA3) {
        xA3 = -yA3 * a + b;
}

void OrangeTriangleSolver::compute_xB3(const double &xA3, double &xB3){
    xB3 = 2 * m_A2B2C2->getABVector().getMiddlePoint().getX() - xA3;
}

void OrangeTriangleSolver::compute_yB3(const double &yA3, double &yB3) {
    yB3 = 2 * m_A2B2C2->getABVector().getMiddlePoint().getY() - yA3;
}

void OrangeTriangleSolver::compute_fg(const double &xA3, const double &yA3, const double &yC3,
                                      double &f, double &g) {
    f = -2 * xA3;
    g = -2 * yA3 * yC3 + pow (yC3, 2)
        - pow(m_A1B1C1->getACVector().getMagnitude(), 2)
        + pow(m_A2B2C2->getACVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getStartPoint().getZ(), 2)
        + pow(xA3, 2)
        + pow(yA3, 2);
}

void OrangeTriangleSolver::compute_hi(const double &xA3, const double &xB3, const double &yA3, const double &yB3,
                                      double &h, double &i) {
    double denominator = -xB3 + xA3;

    i = pow(m_A1B1C1->getBCVector().getMagnitude(), 2)
        - pow(m_A2B2C2->getABVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getEndPoint().getZ(), 2)
        - pow(xB3, 2) - pow(yB3, 2)
        - pow(m_A1B1C1->getACVector().getMagnitude(), 2)
        + pow(m_A2B2C2->getACVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getStartPoint().getZ(), 2)
        + pow(xA3, 2)
        + pow(yA3, 2);

    h = (2 * yB3 - 2 * yA3 ) / (2 * denominator);
    i = i / (2 * denominator);
}

void OrangeTriangleSolver::compute_jkl(const double &xA3, const double &yA3, const double &h, const double &i,
                                       double &j, double &k, double &l) {
    j = pow(h, 2) + 1;
    k = 2 * (h * i - xA3 * h - yA3);
    l = pow(i, 2) - 2 * xA3 * i
        - pow(m_A1B1C1->getACVector().getMagnitude(), 2)
        + pow(m_A2B2C2->getACVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getStartPoint().getZ(), 2)
        + pow(xA3, 2) + pow(yA3, 2);
}

void OrangeTriangleSolver::compute_quadratic_yC3(const double &h, const double &i, const double &j,
                                                 std::vector<double> &yC3) {
    Utils::Math::quadraticEquation(h, i, j, yC3);
}

void OrangeTriangleSolver::compute_xC3(const double &yC3, const double &f, const double &g, double &xC3) {
    xC3 = yC3 * f + g;
}

void OrangeTriangleSolver::compute_yC3(const double &xA3, const double &xB3, const double &yA3, const double &yB3, double &yC3) {
    double q = pow(m_A1B1C1->getBCVector().getMagnitude(), 2)
                - pow(m_A2B2C2->getABVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getEndPoint().getZ(), 2)
                - pow(xB3, 2) - pow(yB3, 2)
                - pow(m_A1B1C1->getACVector().getMagnitude(), 2)
                + pow(m_A2B2C2->getACVector().getEndPoint().getZ() - m_A2B2C2->getACVector().getStartPoint().getZ(), 2)
                + pow(xA3, 2)
                + pow(yA3, 2);

    yC3 = -1 * q / (2 * (yB3 - yA3));
}

void OrangeTriangleSolver::compute_quadratic_xC3(const double &r, const double &f, const double &g,
                                                 std::vector<double> &xC3) {
    Utils::Math::quadraticEquation(1, f, g, xC3);
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