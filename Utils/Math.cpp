//
// Created by Laura Balascuta on 07/03/2021.
//

#include "Math.h"

#include <cmath>

namespace Utils {

void Math::quadraticEquation(const double &a, const double &b, const double &c, std::vector<double> &x) {
    double discriminant = b * b - 4 * a * c;

    if (a != 0 && discriminant > 0) {
        x.push_back((-b + sqrt(discriminant)) / (2 * a));
        x.push_back((-b - sqrt(discriminant)) / (2 * a));
    }

    if (a != 0 && discriminant == 0) {
        x.push_back(-b / (2 * a));
    }
}

}