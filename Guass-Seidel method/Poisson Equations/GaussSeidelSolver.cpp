// GaussSeidelSolver.cpp
#include "GaussSeidelSolver.h"
#include <cmath>
#include <iostream>

GaussSeidelSolver::GaussSeidelSolver(int size) {
    // Constructor (if needed)
}

GaussSeidelSolver::~GaussSeidelSolver() {
    // Destructor (if needed)
}

void GaussSeidelSolver::solve(std::vector<double>& phi, const std::vector<double>& f, double tolerance, int maxIterations) {
    int n = phi.size();
    int iterations = 0;
    double error = tolerance + 1.0;

    while (iterations < maxIterations && error > tolerance) {
        error = 0.0;

        for (int i = 1; i < n - 1; ++i) {
            double oldPhi = phi[i];
            phi[i] = 0.5 * (phi[i - 1] + phi[i + 1] - f[i]);
            error += std::abs(phi[i] - oldPhi);
        }

        ++iterations;
    }

    std::cout << "Gauss-Seidel method converged in " << iterations << " iterations." << std::endl;
}
