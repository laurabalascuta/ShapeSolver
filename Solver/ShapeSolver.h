//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_SHAPESOLVER_H
#define SHAPESOLVER_SHAPESOLVER_H

#include "../Entities/Entity.h"

namespace Solver {

/// \brief Base class for shape problems/solvers
class ShapeSolver {
public:
    /// \brief Default destructor
    virtual ~ShapeSolver() = default;

    /// \brief It solves the given problem
    virtual void solve() = 0;
};

}

#endif //SHAPESOLVER_SHAPESOLVER_H
