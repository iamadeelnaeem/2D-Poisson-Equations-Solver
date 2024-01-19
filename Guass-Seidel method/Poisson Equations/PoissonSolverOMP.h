// PoissonSolverOMP.h
#ifndef PASSION_SOLVER_OMP_H
#define PASSION_SOLVER_OMP_H

#pragma once
#include <vector>

class PoissonSolverOMP {
public:
    PoissonSolverOMP(int size);

    ~PoissonSolverOMP();

    void solve(std::vector<double>& phi, const std::vector<double>& f, double tolerance, int maxIterations);
    double relaxation(const std::vector<double>& phi, const std::vector<double>& f, const std::vector<double>& residual, int i);
    double vectorInnerProduct(const std::vector<double>& a, const std::vector<double>& b);
    void matrixVectorMultiply(const std::vector<double>& phi, std::vector<double>& result);

    // Other member functions...

private:
    int size;

    // Other private member variables...

    void initialize(std::vector<double>& phi, const std::vector<double>& f);

    // Other private member functions...
};
#endif