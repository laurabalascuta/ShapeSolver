//
// Created by Laura Balascuta on 07/03/2021.
//

#ifndef SHAPESOLVER_ORANGETRIANGLESOLVER_H
#define SHAPESOLVER_ORANGETRIANGLESOLVER_H

#include "ShapeSolver.h"

#include "../Entities/Factory.h"
#include "../Entities/Triangle.h"

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
    /// \param [out] a formula parameter named a
    /// \param [out] b formula parameter named b
    void compute_ab(double &a, double &b);

    /// \brief Calculate c, d, e
    /// \param [in] a formula parameter named a
    /// \param [in] b formula parameter named b
    /// \param [out] c formula parameter named c
    /// \param [out] d formula parameter named d
    /// \param [out] e formula parameter named e
    void compute_cde(const double &a, const double &b, double &c, double &d, double &e);

    /// \brief Calculate m, n
    /// \param [in] yA3 the value of the coordinate yA3
    /// \param [out] m formula parameter named m
    /// \param [out] n formula parameter named n
    void compute_mn(const double &yA3, double &m, double &n);

    /// \brief Calculate yA3
    /// \param [out] yA3 the value of coordinate y of point A3
    void compute_yA3(double &yA3);

    /// \brief Calculate yA3
    /// \param [in] c formula parameter named c
    /// \param [in] d formula parameter named d
    /// \param [in] e formula parameter named e
    /// \param [out] yA3 the values of coordinate y of point A3
    void compute_quadratic_yA3(const double &c, const double &d, const double &e,
                               std::vector<double> &yA3);

    /// \brief Calculate xA3
    /// \param [in] yA3 the value of coordinate y of point A3
    /// \param [in] a formula parameter named a
    /// \param [in] b formula parameter named b
    /// \param [out] xA3 the value of coordinate x of point A3
    void compute_xA3(const double &yA3, const double &a, const double &b, double &xA3);

    /// \brief Calculate xA3
    /// \param [in] m formula parameter named m
    /// \param [in] n formula parameter named n
    /// \param [out] xA3 the values of coordinate x of point A3
    void compute_quadratic_xA3(const double &m, const double &n, std::vector<double> &xA3);

    /// \brief Calculate xB3
    /// \param [in] xA3 the first value of coordinate x of point A3
    /// \param [out] xB3 the value of coordinate x of point B3
    void compute_xB3(const double &xA3, double &xB3);

    /// \brief Calculate yB3
    /// \param [in] yA3 the value of coordinate y of point A3
    /// \param [out] yB3 the value of coordinate y of point B3
    void compute_yB3(const double &yA3, double &yB3);

    /// \brief Calculate f and g
    /// \param [in] xA3 the value of coordinate x of point A3
    /// \param [in] yA3 the value of coordinate y of point A3
    /// \param [in] yC3 the value of coordinate y of point C3
    /// \param [out] f formula parameter named f
    /// \param [out] g formula parameter named g
    void compute_fg(const double &xA3, const double &yA3, const double &yC3,
                    double &f, double &g);

    /// \brief Calculate h and i
    /// \param [in] xA3 the value of coordinate x of point A3
    /// \param [in] xB3 the value of coordinate x of point B3
    /// \param [in] yA3 the value of coordinate y of point A3
    /// \param [in] yB3 value of coordinate y of point B3
    /// \param [out] h formula parameter named h
    /// \param [out] i formula parameter named i
    void compute_hi(const double &xA3, const double &xB3, const double &yA3, const double &yB3,
                    double &h, double &i);

    /// \brief Calculate j, k and l
    /// \param [in] xA3 the value of coordinate x of point A3
    /// \param [in] yA3 the value of coordinate y of point A3
    /// \param [in] h formula parameter named h
    /// \param [in] i formula parameter named i
    /// \param [out] j formula parameter named j
    /// \param [out] k formula parameter named k
    /// \param [out] l formula parameter named l
    void compute_jkl(const double &xA3, const double &yA3, const double &h, const double &i,
                     double &j, double &k, double &l);


    /// \brief Calculate yC3
    /// \param [in] xA3 the value of coordinate x of point A3
    /// \param [in] xB3 the value of coordinate x of point B3
    /// \param [in] yA3 the value of coordinate y of point A3
    /// \param [in] yB3 the value of coordinate y of point B3
    /// \param [out] yC3 the value of coordinate y of point C3
    void compute_yC3(const double &xA3, const double &xB3, const double &yA3, const double &yB3, double &yC3);

    /// \brief Calculate yC3
    /// \param [in] h formula parameter named h
    /// \param [in] i formula parameter named i
    /// \param [in] j formula parameter named j
    /// \param [out] yC3 the values of coordinate y of point C3
    void compute_quadratic_yC3(const double &h, const double &i, const double &j,
                               std::vector<double> &yC3);

    /// \brief Calculate xC3
    /// \param [in] yC3 the value of coordinate y of point C3
    /// \param [in] f formula parameter named f
    /// \param [in] g formula parameter named g
    /// \param [out] xC3 the value of coordinate x of point C3
    void compute_xC3(const double &yC3, const double &f, const double &g, double &xC3);

    /// \brief Calculate xC3
    /// \param [in] r formula parameter named r
    /// \param [in] f formula parameter named f
    /// \param [in] g formula parameter named g
    /// \param [out] xC3 the value of coordinate x of point C3
    void compute_quadratic_xC3(const double &r, const double &f, const double &g,
                               std::vector<double> &xC3);

    /// \brief Call the assumptions check assumptions, the creation of the triangle and add it to the list
    /// \param [in] coord the coordinates of the triangle
    /// \param [in] name the name of the triangle
    void addTriangleToOutputList(const std::vector<double> &coord, const std::vector<std::string> &name);

    /// \brief Check all the assumptions of the problem
    [[nodiscard]] bool checkAssumptionsSolved(const std::shared_ptr<Entities::Triangle> &triangle) const;

    /// \brief Check assumption 1: A2B2C2 and A3B3C3 are in the same plane
    /// \return true if checked
    [[nodiscard]] bool checkSamePlane(const std::shared_ptr<Entities::Triangle> &triangle) const;

    /// \brief Check assumption 2: A1B1C1 and A3B3C3 have the same vector magnitudes
    /// \return true if checked
    [[nodiscard]] bool checkMagnitudeEquality(const std::shared_ptr<Entities::Triangle> &triangle) const;

    /// \brief Check assumption 3: A1B1 is perpendicular on A3B3
    /// \return true if checked
    [[nodiscard]] bool checkPerpendicularVectors(const std::shared_ptr<Entities::Triangle> &triangle) const;

    /// \brief Check assumption 4: middle of A2B2 is the same as the middle of A3B3
    /// \return true if checked
    [[nodiscard]] bool checkMiddle(const std::shared_ptr<Entities::Triangle> &triangle) const;

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

    bool m_solved;
};

}

#endif //SHAPESOLVER_ORANGETRIANGLESOLVER_H
