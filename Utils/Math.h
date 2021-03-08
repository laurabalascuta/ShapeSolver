//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_MATH_H
#define SHAPESOLVER_MATH_H

#include <optional>

namespace Utils {

// Class that holds mathematical utilities
class Math {
public:
    /// \brief Solves a quadratic equation
    /// \param [in] a the quadratic coefficient
    /// \param [in] b the linear coefficient
    /// \param [in] c the constant
    /// \param [out] x1 the first solution of the equation
    /// \param [out] x2 the second solution of the equation
    static void quadraticEquation(const double &a, const double &b, const double &c,
                                  std::optional<double> &x1, std::optional<double> &x2);
};

}

#endif //SHAPESOLVER_MATH_H
