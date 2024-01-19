// MultigridSolver.h
#ifndef MULTIGRID_SOLVER_H
#define MULTIGRID_SOLVER_H

#include <vector>

class MultigridSolver {
public:
    MultigridSolver(int size);
    ~MultigridSolver();
    void solve(std::vector<double>& phi, const std::vector<double>& f, double tolerance, int maxIterations);
    void restriction(const std::vector<double>& r, std::vector<double>& r2h, int n);
    void interpolation(const std::vector<double>& E2h, std::vector<double>& E, int n);
    std::vector<double> cycle(const std::vector<double>& r, int lambda, int theta);
    std::vector<double> gaussSeidel(const std::vector<double>& x, const std::vector<double>& b);
    std::vector<double> computeResidual(const std::vector<double>& b, const std::vector<double>& x);
    double computeNorm(const std::vector<double>& vec);
};

#endif // MULTIGRID_SOLVER_H
