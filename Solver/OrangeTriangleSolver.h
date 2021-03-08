//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_ORANGETRIANGLESOLVER_H
#define SHAPESOLVER_ORANGETRIANGLESOLVER_H

#include "ShapeSolver.h"

#include "../Entities/Factory.h"
#include "../Entities/Triangle.h"

#include <optional>
#include <vector>

namespace Solver {

/// \brief Class for the orange triangle problem
///
/// The problem states the following:
/// <br> - input: A1B1C1 green triangle and A2B2C2 blue triangle
/// <br> - output: A3B3C3 orange triangle
/// <br> - assumptions:
/// <br>  1. A2B2C2 and A3B3C3 are in the same plane
/// <br>  2. A1B1C1 and A3B3C3 have the same vector magnitudes
/// <br>  3. A1B1 is perpendicular on A3B3
/// <br>  4. middle of A2B2 is the same as the middle of A3B3

class OrangeTriangleSolver : public  ShapeSolver {
public:
    /// \brief Default constructor
    OrangeTriangleSolver();

    /// \brief Default destructor
    ~OrangeTriangleSolver() override = default;

    /// \brief It solves the orange triangle problem
    void solve() override;

private:
    /// \brief Calculate a and b
    /// \param [out] a
    /// \param [out] b
    void ab(double &a, double &b);

    /// \brief Calculate c, d, e
    /// \param [in] a
    /// \param [in] b
    /// \param [out] c
    /// \param [out] d
    /// \param [out] e
    void cde(const double &a, const double &b, double &c, double &d, double &e);

    /// \brief Calculate yA3
    /// \param [in] c
    /// \param [in] d
    /// \param [in] e
    /// \param [out] yA3_1 the first possible value of yA3
    /// \param [out] yA3_2 the second possible value of yA3
    void yA3(const double &c, const double &d, const double &e,
             std::optional<double> &yA3_1,std::optional<double> &yA3_2);

    /// \brief Calculate xA3
    /// \param [in] yA3_1
    /// \param [in] yA3_2
    /// \param [in] a
    /// \param [in] b
    /// \param [out] xA3_1 the first possible value of xA3
    /// \param [out] xA3_2 the first possible value of xA3
    void xA3(const std::optional<double> &yA3_1, const std::optional<double> &yA3_2,
             const double &a, const double &b,
             std::optional<double> &xA3_1, std::optional<double> &xA3_2);

    /// \brief Calculate xB3
    /// \param [in] xA3_1
    /// \param [in] xA3_2
    /// \param [out] xB3_1 the first possible value of xB3
    /// \param [out] xB3_2 the second possible value of xB3
    void xB3(const std::optional<double> &xA3_1, const std::optional<double> &xA3_2,
             std::optional<double> &xB3_1, std::optional<double> &xB3_2);

    /// \brief Calculate yB3
    /// \param [in] yA3_1
    /// \param [in] yA3_2
    /// \param [out] yB3_1 the first possible value of yB3
    /// \param [out] yB3_2 the second possible value of yB3
    void yB3(const std::optional<double> &yA3_1, const std::optional<double> &yA3_2,
             std::optional<double> &yB3_1, std::optional<double> &yB3_2);

    /// \brief Calculate f and g
    /// \param [in] xA3_1
    /// \param [in] xB3_1
    /// \param [in] yA3_1
    /// \param [in] yB3_1
    /// \param [out] f
    /// \param [out] g
    void fg(const std::optional<double> &xA3_1, const std::optional<double> &xB3_1,
            const std::optional<double> &yA3_1, const std::optional<double> &yB3_1,
            double &f, double &g);

    /// \brief Calculate h, i and j
    /// \param [in] xA3_1
    /// \param [in] yA3_1
    /// \param [in] f
    /// \param [in] g
    /// \param [out] h
    /// \param [out] i
    /// \param [out] j
    void hij(const std::optional<double> &xA3_1, const std::optional<double> &yA3_1,
             const double &f, const double &g,
             double &h, double &i, double &j);

    /// \brief Calculate yC3
    /// \param [in] h
    /// \param [in] i
    /// \param [in] j
    /// \param [out] yC3_1 the first possible value of yC3
    /// \param [out] yC3_2 the second possible value of yC3
    void yC3(const double &h, const double &i, const double &j,
             std::optional<double> &yC3_1, std::optional<double> &yC3_2);

    /// \brief Calculate xC3
    /// \param [in] yC3_1
    /// \param [in] yC3_2
    /// \param [in] f
    /// \param [in] g
    /// \param [out] xC3_1 the first possible value of xC3
    /// \param [out] xC3_2 the second possible value of xC3
    void xC3(const std::optional<double> &yC3_1, const std::optional<double> &yC3_2,
             const double &f, const double &g,
             std::optional<double> &xC3_1, std::optional<double> &xC3_2);

    /// \brief Print the input triangles
    void printInput() const;

    /// \brief Print the output triangles
    void printOutput() const;

private:
    static const int NO_INPUT_TRIANGLES;
    static const int NO_COORD_PER_POINT;
    static const std::string COORDINATES_FILENAME;

    std::shared_ptr<Entities::Triangle> m_A1B1C1; // input - green triangle
    std::shared_ptr<Entities::Triangle> m_A2B2C2; // input - blue triangle

    std::vector<std::shared_ptr<Entities::Triangle>> m_outputTriangles;
    std::unique_ptr<Entities::Factory> m_factory;

    bool m_solved = false;
};

}

#endif //SHAPESOLVER_ORANGETRIANGLESOLVER_H
