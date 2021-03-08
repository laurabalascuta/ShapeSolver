//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_POINT_H
#define SHAPESOLVER_POINT_H

#include "Entity.h"

#include <iostream>

namespace Entities {

/// \brief Class that models a point
class Point : public Entity {
public:
    /// \brief Default constructor is deleted to prevent empty creation
    Point() = delete;

    /// \brief Parametrized constructor
    /// \param [in] x the coordinate x
    /// \param [in] y the coordinate y
    /// \param [in] z the coordinate y
    /// \param [in] name the name of the point
    Point(const double &x, const double &y, const double &z, const std::string &name);

    /// \brief Default destructor
    ~Point() override = default;

    /// Returns the x coordinate of the point
    /// \return the x coordinate
    [[nodiscard]] double getX() const;

    /// Returns the y coordinate of the point
    /// \return the y coordinate
    [[nodiscard]] double getY() const;

    /// Returns the z coordinate of the point
    /// \return the z coordinate
    [[nodiscard]] double getZ() const;

    /// \brief It prints the coordinates and the name that defines the point
    ///        using the output stream insertion operator
    void printInfo() const override;

    /// \brief It returns the name of the point
    /// \return the name of the point
    [[nodiscard]] std::string getName() const override;

    /// Overload of the output stream insertion operator
    /// \param os the ostream object to write data to
    /// \param point the point object to write out
    /// \return the ostream object
    friend std::ostream &operator<<(std::ostream &os, const Point &point);

private:
    double m_x;
    double m_y;
    double m_z;

    std::string m_name;
};

}

#endif //SHAPESOLVER_POINT_H
