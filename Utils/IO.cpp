//
// Created by Laura Balascuta on 07/03/2021.
//

#include "IO.h"

#include <iostream>
#include <fstream>
#include <string>

namespace Utils {

std::vector<double> IO::readCoordinatesFromFile(const std::string& fileName) {
    std::vector<double> coordinates;

    std::ifstream file(fileName);

    double x, y, z;
    while (file >> x >> y >> z)
    {
        coordinates.push_back(x);
        coordinates.push_back(y);
        coordinates.push_back(z);
    }

    return coordinates;
}

}