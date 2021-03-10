#include "Solver/OrangeTriangleSolver.h"
#include "Utils/Math.h"

int main() {

  std::unique_ptr<Solver::ShapeSolver> solver(new Solver::OrangeTriangleSolver());
  solver->solve();

  return 0;
}
