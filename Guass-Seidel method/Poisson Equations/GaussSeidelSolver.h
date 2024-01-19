// GaussSeidelSolver.h
#ifndef GAUSS_SEIDEL_SOLVER_H
#define GAUSS_SEIDEL_SOLVER_H

#include <vector>

class GaussSeidelSolver {
public:
    GaussSeidelSolver(int size);
    ~GaussSeidelSolver();
    void solve(std::vector<double>& phi, const std::vector<double>& f, double tolerance, int maxIterations);
};

#endif // GAUSS_SEIDEL_SOLVER_H
