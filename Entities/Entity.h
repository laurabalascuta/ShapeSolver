//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_ENTITY_H
#define SHAPESOLVER_ENTITY_H

#include <ostream>

namespace Entities {

/// \brief The base class for an Entity
class Entity {
public:

    /// \brief Default destructor
    virtual ~Entity() = default;

    /// \brief It prints the information that defines the Entity
    virtual void printInfo() const = 0;

    /// \brief It returns the name of an Entity
    /// \return the name of an Entity
    virtual std::string getName() const = 0;
};

}

#endif //SHAPESOLVER_ENTITY_H
