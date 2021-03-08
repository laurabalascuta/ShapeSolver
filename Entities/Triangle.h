//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_TRIANGLE_H
#define SHAPESOLVER_TRIANGLE_H

#include "Entity.h"
#include "Vector.h"

namespace Entities {

class Triangle : public Entity {
public:
    static const int POINTS_NO;

public:
    /// \brief Default constructor is deleted to prevent empty creation
    Triangle() = delete;

    /// \brief Parametrized constructor
    /// \param [in] AB the vector AB
    /// \param [in] BC the vector BC
    /// \param [in] AC the vector AC
    Triangle(const Vector& AB, const Vector& BC, const Vector& AC);

    /// \brief Default destructor
    ~Triangle();

    /// \brief Get the AB vector of the triangle
    /// \return the AB vector of the triangle
    Vector getABVector() const;

    /// \brief Get the BC vector of the triangle
    /// \return the BC vector of the triangle
    Vector getBCVector() const;

    /// \brief Get the AC vector of the triangle
    /// \return the AC vector of the triangle
    Vector getACVector() const;

    /// \brief It prints the coordinates of the points and the name of the Triangle
    ///        using the output stream insertion operator
    void printInfo() const override;

    /// \brief It returns the name of the Triangle
    /// \return the name of the Triangle
    std::string getName() const override;

    /// Overload of the output stream insertion operator
    /// \param os the ostream object to write data to
    /// \param point the Triangle object to write out
    /// \return the ostream object
    friend std::ostream &operator<<(std::ostream &os, const Triangle &triangle);

private:
    Vector m_AB;
    Vector m_BC;
    Vector m_AC;

    std::string m_name;
};

}

#endif //SHAPESOLVER_TRIANGLE_H
