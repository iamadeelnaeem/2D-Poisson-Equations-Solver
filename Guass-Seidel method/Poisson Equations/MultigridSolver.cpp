#include "MultigridSolver.h"
#include <cmath>
#include <iostream>

MultigridSolver::MultigridSolver(int size) {
    // Constructor (if needed)
}

MultigridSolver::~MultigridSolver() {
    // Destructor (if needed)
}

void MultigridSolver::solve(std::vector<double>& phi, const std::vector<double>& f, double tolerance, int maxIterations) {
    int n = std::sqrt(phi.size()) - 1; // Assuming n is a power of 2
    int lambda = 1;
    int theta = 1;
    int dim = phi.size();

    // Set up matrices (A) - you need to define this function
    // matrixSetup(dim);

    // Perform initial relaxation
    for (int i = 0; i < lambda; ++i) {
        // Relaxation using Gauss-Seidel or any other smoother
        phi = gaussSeidel(phi, f);
    }

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        std::vector<double> residual(dim, 0.0);

        // Compute the residual: residual = f - A * phi
        residual = computeResidual(f, phi);

        // Compute the norm of the residual: residualNorm = sqrt(residual * residual)
        double residualNorm = computeNorm(residual);

        if (residualNorm < tolerance) {
            std::cout << "Multigrid method converged in " << iteration << " iterations." << std::endl;
            return;
        }

        // Perform a Multigrid cycle
        phi = cycle(residual, lambda, theta);
    }

    std::cout << "Multigrid method did not converge within the specified number of iterations." << std::endl;
}

double MultigridSolver::computeNorm(const std::vector<double>& vec) {
    double norm = 0.0;

    for (double val : vec) {
        norm += val * val;
    }

    return std::sqrt(norm);
}

void MultigridSolver::restriction(const std::vector<double>& r, std::vector<double>& r2h, int n) {
    // Implementation of restriction
    int N2h = (n + 1) / 2 - 1;
    int dim2h = N2h * N2h;
    r2h.assign(dim2h, 0.0);

    for (int i = 1, l = 0, k = 0; i <= n; ++i) {
        for (int j = 1; j <= n; ++j, ++k) {
            if (i % 2 == 0 && j % 2 == 0) {
                r2h[l] = 1.0 / 16.0 * (4.0 * r[k] + 2.0 * (r[k - 1] + r[k + 1] + r[k - n] + r[k + n])
                    + r[k + n - 1] + r[k + n + 1] + r[k - n - 1] + r[k - n + 1]);
                ++l;
            }
        }
    }
}

void MultigridSolver::interpolation(const std::vector<double>& E2h, std::vector<double>& E, int n) {
    // Implementation of interpolation
    int dim = E.size();
    E.assign(dim, 0.0);

    for (int i = 1, k = 0, l = 0; i <= n; ++i) {
        for (int j = 1; j <= n; ++j, ++k) {
            if (i % 2 == 0 && j % 2 == 0) {
                E[k] = E2h[l];
                ++l;
            }
        }
    }

    for (int i = 1, k = 0; i <= n; ++i) {
        for (int j = 1; j <= n; ++j, ++k) {
            if (i % 2 == 0 && j % 2 != 0) {
                if (j != 1 && j != n) E[k] = 1.0 / 2.0 * (E[k - 1] + E[k + 1]);
                if (j == 1) E[k] = 1.0 / 1.0 * (E[k + 1]);
                if (j == n) E[k] = 1.0 / 1.0 * (E[k - 1]);
            }
            if (i % 2 != 0 && j % 2 == 0) {
                if (i != 1 && i != n) E[k] = 1.0 / 2.0 * (E[k - n] + E[k + n]);
                if (i == 1) E[k] = 1.0 / 1.0 * (E[k + n]);
                if (i == n) E[k] = 1.0 / 1.0 * (E[k - n]);
            }
            if (i % 2 != 0 && j % 2 != 0) {
                if (i != 1 && j != 1 && i != n && j != n) E[k] = 1.0 / 4.0 * (E[k + n - 1] + E[k + n + 1] + E[k - n + 1] + E[k - n - 1]);
                if (i == 1 && j == 1) E[k] = 1.0 / 1.0 * (E[k + n + 1]);
                if (i == n && j == n) E[k] = 1.0 / 1.0 * (E[k - n - 1]);
                if (i == 1 && j == n) E[k] = 1.0 / 1.0 * (E[k + n - 1]);
                if (i == n && j == 1) E[k] = 1.0 / 1.0 * (E[k - n + 1]);
                if (i == 1 && j != 1 && j != n) E[k] = 1.0 / 2.0 * (E[k + n + 1] + E[k + n - 1]);
                if (i == n && j != 1 && j != n) E[k] = 1.0 / 2.0 * (E[k - n - 1] + E[k - n + 1]);
                if (j == 1 && i != 1 && i != n) E[k] = 1.0 / 2.0 * (E[k + n + 1] + E[k - n + 1]);
                if (j == n && i != 1 && i != n) E[k] = 1.0 / 2.0 * (E[k + n - 1] + E[k - n - 1]);
            }
        }
    }
}

std::vector<double> MultigridSolver::cycle(const std::vector<double>& r, int lambda, int theta) {
    int dim = r.size();
    int n = std::sqrt(dim) - 1; // Assuming n is a power of 2

    std::vector<double> x(dim, 0.0);
    std::vector<double> residual(r);

    // Pre-smoothing
    for (int i = 0; i < lambda; ++i) {
        // Relaxation using Gauss-Seidel or any other smoother
        x = gaussSeidel(x, residual);
    }

    // Compute residual
    // residual = r - A * x
    // You need to define matrix-vector multiplication (A * x)
    residual = computeResidual(r, x);

    // Restriction
    std::vector<double> r2h;
    restriction(residual, r2h, n);

    // Solve error equation A2h * E2h = r2h recursively using the same Multigrid solver
    int N2h = (n + 1) / 2 - 1;
    int dim2h = N2h * N2h;
    std::vector<double> E2h(dim2h, 0.0);
    E2h = cycle(r2h, lambda, theta);

    // Interpolation
    std::vector<double> E(dim, 0.0);
    interpolation(E2h, E, n);

    // Update the solution
    for (int i = 0; i < dim; ++i) {
        x[i] += E[i];
    }

    // Post-smoothing
    for (int i = 0; i < theta; ++i) {
        // Relaxation using Gauss-Seidel or any other smoother
        x = gaussSeidel(x, residual);
    }

    return x;
}

std::vector<double> MultigridSolver::gaussSeidel(const std::vector<double>& x, const std::vector<double>& b) {
    int n = std::sqrt(x.size()) - 1;
    std::vector<double> result(x);

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            int index = i * (n + 1) + j;
            result[index] = 0.25 * (result[index - 1] + result[index + 1]
                + result[index - (n + 1)] + result[index + (n + 1)]
                - b[index]);
        }
    }

    return result;
}

std::vector<double> MultigridSolver::computeResidual(const std::vector<double>& b, const std::vector<double>& x) {
    int n = std::sqrt(b.size()) - 1;
    std::vector<double> residual(b);

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            int index = i * (n + 1) + j;
            residual[index] -= (x[index - 1] + x[index + 1]
                + x[index - (n + 1)] + x[index + (n + 1)]);
        }
    }

    return residual;
}
