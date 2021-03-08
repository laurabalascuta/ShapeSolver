//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_MATH_H
#define SHAPESOLVER_MATH_H

#include <vector>

namespace Utils {

// Class that holds mathematical utilities
class Math {
public:
    /// \brief Solves a quadratic equation
    /// \param [in] a the quadratic coefficient
    /// \param [in] b the linear coefficient
    /// \param [in] c the constant
    /// \param [out] x the roots
    static void quadraticEquation(const double &a, const double &b, const double &c, std::vector<double> &x);
};

}

#endif //SHAPESOLVER_MATH_H
