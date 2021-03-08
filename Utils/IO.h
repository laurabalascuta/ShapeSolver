//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_IO_H
#define SHAPESOLVER_IO_H

#include <string>
#include <vector>

namespace Utils {

// \brief Class that holds IO utilities
class IO {
public:
    /// \brief Reads the coordinates from a file and returns them. <br>
    ///        The x, y, z coordinates of a point must lay on a single line separated by spaces.<br>
    ///         1 2 3 <br>
    ///         4 5 6
    /// \param [in] fileName the name of the file to be read
    /// \return the coordinates
    static std::vector<double> readCoordinatesFromFile(const std::string& fileName);
};

}

#endif //SHAPESOLVER_IO_H
