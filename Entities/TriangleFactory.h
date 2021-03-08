//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_TRIANGLEFACTORY_H
#define SHAPESOLVER_TRIANGLEFACTORY_H

#include "Factory.h"
#include "Triangle.h"

namespace Entities {

/// \brief Class that creates triangles
class TriangleFactory : public Factory {
public:
    /// \brief Default destructor
    virtual ~TriangleFactory() = default;

    /// \brief Creates an Entity from the input data
    /// \param [in] coord the coordinates of the Entity
    /// \param [in] name the name of the Entity
    /// \return the new Entity
    std::shared_ptr<Entity> createEntity(const std::vector<double> &coord, const std::vector<std::string> name) const override;
};

}

#endif //SHAPESOLVER_TRIANGLEFACTORY_H
