// main.cpp
#include "GaussSeidelSolver.h"
#include <iostream>
#include <vector>

int main() {
    int n;
    std::cout << "Enter the number of grid points: ";
    std::cin >> n;

    std::vector<double> phi(n, 0.0);
    std::vector<double> f(n, 0.0);

    std::cout << "Enter coefficients for the Poisson equation:" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << "Enter f[" << i << "]: ";
        std::cin >> f[i];
    }

    double tolerance = 1e-6;
    int maxIterations = 10000;

    // Create an instance of GaussSeidelSolver
    GaussSeidelSolver gaussSeidelSolver(n);

    // Solve the Poisson equation using Gauss-Seidel method
    gaussSeidelSolver.solve(phi, f, tolerance, maxIterations);

    std::cout << "Final solution (phi): ";
    for (int i = 0; i < n; ++i) {
        std::cout << phi[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
