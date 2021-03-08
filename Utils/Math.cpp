//
// Created by Laura Balascuta on 07/03/2021.
//

#include "Math.h"

#include <math.h>

namespace Utils {

void Math::quadraticEquation(const double &a, const double &b, const double &c,
                             std::optional<double> &x1, std::optional<double> &x2) {
    double discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {
        x1 = (-b + sqrt(discriminant)) / (2 * a);
        x2 = (-b - sqrt(discriminant)) / (2 * a);
    }

    if (discriminant == 0) {
        x1 = -b / (2 * a);
    }
}

}