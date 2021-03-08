//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_VECTOR_H
#define SHAPESOLVER_VECTOR_H

#include "Entity.h"
#include "Point.h"

#include <iostream>

namespace Entities {

class Vector : public Entity {
public:
    /// \brief Default constructor is deleted to prevent empty creation
    Vector() = delete;

    /// \brief Parametrized constructor
    /// \param [in] start the start point of the vector
    /// \param [in] end the end point of the vector
    Vector(const Point &start, const Point &end);

    /// \brief Default destructor
    ~Vector() override = default;

    /// \brief Get the start point of the vector
    /// \return the start point of the vector
    [[nodiscard]] Point getStartPoint() const;

    /// \brief Get the end point of the vector
    /// \return the end point of the vector
    [[nodiscard]] Point getEndPoint() const;

    /// \brief Get the position vector's end point of the current vector
    ///        The position vector is the position of a point P in space
    ///        in relation to an arbitrary reference origin O.
    ///        It is obtained by subtracting the coordinates of the initial point
    ///        from the coordinates of the terminal point.
    ///        Formula: Point(x_end - x_start, y_end - y_start, z_end - z_start)
    /// \return the position vector's end point
    [[nodiscard]] Point getPositionVectorPoint() const;

    /// \brief Get the middle point of the vector
    ///        The middle point has as coordinates the start plus end divided by 2.
    ///        Formula: Point((x_start + x_end)/2, (y_start + y_end)/2, (z_start + z_end)/2)
    /// \return the middle point of the vector
    [[nodiscard]] Point getMiddlePoint() const;

    /// \brief Get the magnitude of the vector
    ///        Formula: sqrt(pow((x_end - x_start), 2) + pow((y_end - y_start), 2) + pow((z_start - z_end), 2));
    /// \return the magnitude of the vector
    [[nodiscard]] double getMagnitude() const;

    /// \brief Compute the dot product of this vector and an input vector
    /// \param [in] vector the input vector
    /// \return the value of the dot product
    [[nodiscard]] double computeDotProduct(const Vector& vector) const;

    /// \brief It prints the start and end points coordinates and the name of the Vector
    ///        using the output stream insertion operator
    void printInfo() const override;

    /// \brief It returns the name of the vector
    /// \return the name of the vector
    [[nodiscard]] std::string getName() const override;

    /// Overload of the output stream insertion operator
    /// \param os the ostream object to write data to
    /// \param point the vector object to write out
    /// \return the ostream object
    friend std::ostream &operator<<(std::ostream &os, const Vector &vector);

private:
    static const std::string MIDDLE_POINT_NAME;

    Point m_start;
    Point m_end;

    Point m_positionVector;
    Point m_middle;

    std::string m_name;

    double m_magnitude;
};

}

#endif //SHAPESOLVER_VECTOR_H
