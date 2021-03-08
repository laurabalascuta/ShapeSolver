//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_FACTORY_H
#define SHAPESOLVER_FACTORY_H

#include "Entity.h"

#include <memory>
#include <vector>

namespace Entities {

/// \brief The base class for the Entity Factory
class Factory {
public:
    /// \brief Default destructor
    virtual ~Factory() = default;

    /// \brief Creates an Entity from the input data
    /// \param [in] coord the coordinates of the Entity
    /// \param [in] name the name of the Entity
    /// \return the new Entity
    virtual std::shared_ptr<Entity> createEntity(const std::vector<double> &coord, const std::vector<std::string> name) const = 0;
};

}
#endif //SHAPESOLVER_FACTORY_H
