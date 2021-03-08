//
// Created by Laura Balascuta on 07/03/2021.
//

#include "TriangleFactory.h"

namespace Entities {

std::shared_ptr<Entity> TriangleFactory::createEntity(const std::vector<double> &coord, const std::vector<std::string> &name) const {
    std::shared_ptr<Triangle> triangle = nullptr;

    // check if there are enough coordinates and names in the input to create a vector
    if (coord.size() == Triangle::POINTS_NO * 3) {
        if (name.size() == Triangle::POINTS_NO) {
            Point A1(coord[0], coord[1], coord[2], name[0]);
            Point B1(coord[3], coord[4], coord[5], name[1]);
            Point C1(coord[6], coord[7], coord[8], name[2]);

            Vector A1B1(A1, B1);
            Vector B1C1(B1, C1);
            Vector A1C1(A1, C1);

            // check if the result would be a valid triangle
            bool valid = A1B1.getMagnitude() + B1C1.getMagnitude() > A1C1.getMagnitude() &&
                B1C1.getMagnitude() + A1C1.getMagnitude() > A1B1.getMagnitude() &&
                A1B1.getMagnitude() + A1C1.getMagnitude() > B1C1.getMagnitude();
            if (valid) {
                triangle = std::shared_ptr<Triangle>(new Triangle(A1B1, B1C1, A1C1));
            } else {
                std::cout << "Factory: The provided coordinates for triangle "
                          << name[0] + name[1] + name[2]
                          << " are not valid to form a triangle." << std::endl;
            }
        }

    } else {
        std::cout << "Factory: There are not enough coordinates provided to create a triangle." << std::endl;
    }

    return triangle;
}

}