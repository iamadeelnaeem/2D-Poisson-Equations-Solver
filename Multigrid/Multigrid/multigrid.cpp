// multigrid.cpp

#include "multigrid.h"
#include <iostream>

MultigridSolver::MultigridSolver(int size) : N(size) {
    // Initialize the grids and other member variables
    u = new double* [N];
    rhs = new double* [N];
    for (int i = 0; i < N; ++i) {
        u[i] = new double[N]();
        rhs[i] = new double[N]();
    }
    // Additional initialization if needed
}

MultigridSolver::~MultigridSolver() {
    // Deallocate memory
    for (int i = 0; i < N; ++i) {
        delete[] u[i];
        delete[] rhs[i];
    }
    delete[] u;
    delete[] rhs;
    // Additional cleanup if needed
}

void MultigridSolver::printGrid() const {
    // Print the resulting grid
    std::cout << "Resulting Grid:" << std::endl;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << u[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Printed grid successfully." << std::endl;
}


void MultigridSolver::setGrid(double** initialGuess) {
    // Set the initial guess for the grid
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] = initialGuess[i][j];
        }
    }
}

double MultigridSolver::getSolution(int i, int j) {
    // Return the solution at position (i, j)
    return u[i][j];
}

void MultigridSolver::initializeRHS() {
    // Implementation to initialize the right-hand side
    // For simplicity, setting RHS to a constant value
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            rhs[i][j] = 1.0;
        }
    }
}

void MultigridSolver::initializeGrid() {
    // Initialize the grid with zeros
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u[i][j] = 0.0;
        }
    }
}

void MultigridSolver::setBoundaryValues() {
    // Ask the user to input grid points on the boundary
    std::cout << "Enter values for grid points on the boundary:" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Boundary point (" << i << ", 0): ";
        std::cin >> u[i][0];
        std::cout << "Boundary point (" << i << ", " << N - 1 << "): ";
        std::cin >> u[i][N - 1];
    }
    for (int j = 0; j < N; ++j) {
        std::cout << "Boundary point (0, " << j << "): ";
        std::cin >> u[0][j];
        std::cout << "Boundary point (" << N - 1 << ", " << j << "): ";
        std::cin >> u[N - 1][j];
    }
}

void MultigridSolver::smooth() {
    // Implementation of a basic Jacobi smoother
    double** temp = new double* [N];
    for (int i = 0; i < N; ++i) {
        temp[i] = new double[N]();
    }

    for (int iter = 0; iter < 5; ++iter) {
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                temp[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - rhs[i][j]);
            }
        }

        // Update 'u' with the temporary values
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                u[i][j] = temp[i][j];
            }
        }
    }

    // Deallocate temporary array
    for (int i = 0; i < N; ++i) {
        delete[] temp[i];
    }
    delete[] temp;
}

void MultigridSolver::restrict() {
    // Implementation of a basic restriction operator (simple averaging)
    int coarseSize = N / 2;
    double** coarse = new double* [coarseSize];
    for (int i = 0; i < coarseSize; ++i) {
        coarse[i] = new double[coarseSize]();
    }

    for (int i = 0; i < coarseSize; ++i) {
        for (int j = 0; j < coarseSize; ++j) {
            coarse[i][j] = 0.25 * (u[2 * i][2 * j] + u[2 * i + 1][2 * j] + u[2 * i][2 * j + 1] + u[2 * i + 1][2 * j + 1]);
        }
    }

    // Replace 'u' with the coarser grid values
    for (int i = 0; i < coarseSize; ++i) {
        for (int j = 0; j < coarseSize; ++j) {
            u[i][j] = coarse[i][j];
        }
    }

    // Deallocate memory for the coarse grid
    for (int i = 0; i < coarseSize; ++i) {
        delete[] coarse[i];
    }
    delete[] coarse;
}

void MultigridSolver::prolongate() {
    // Implementation of a basic prolongation operator (linear interpolation)
    int coarseSize = N / 2;
    double** coarse = new double* [coarseSize];
    for (int i = 0; i < coarseSize; ++i) {
        coarse[i] = new double[coarseSize]();
    }

    // Copy the coarse grid values to the fine grid locations
    for (int i = 0; i < coarseSize; ++i) {
        for (int j = 0; j < coarseSize; ++j) {
            u[2 * i][2 * j] = coarse[i][j];
            u[2 * i + 1][2 * j] = coarse[i][j];
            u[2 * i][2 * j + 1] = coarse[i][j];
            u[2 * i + 1][2 * j + 1] = coarse[i][j];
        }
    }

    // Deallocate memory for the coarse grid
    for (int i = 0; i < coarseSize; ++i) {
        delete[] coarse[i];
    }
    delete[] coarse;
}

void MultigridSolver::solve() {
    // Number of multigrid cycles (adjust as needed)
    int numCycles = 5;

    for (int cycle = 0; cycle < numCycles; ++cycle) {
        // Pre-smoothing
        smooth();

        // Compute the residual
        double** residual = new double* [N];
        for (int i = 0; i < N; ++i) {
            residual[i] = new double[N]();
        }
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                residual[i][j] = rhs[i][j] - (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]);
            }
        }

        // Print the residual grid (add this part)
        std::cout << "Residual Grid (Cycle " << cycle + 1 << "):" << std::endl;
        printGrid();

        // Restriction
        int coarseSize = N / 2;
        double** coarseResidual = new double* [coarseSize];
        for (int i = 0; i < coarseSize; ++i) {
            coarseResidual[i] = new double[coarseSize]();
        }
        for (int i = 0; i < coarseSize; ++i) {
            for (int j = 0; j < coarseSize; ++j) {
                coarseResidual[i][j] = 0.25 * (residual[2 * i][2 * j] + residual[2 * i + 1][2 * j] +
                    residual[2 * i][2 * j + 1] + residual[2 * i + 1][2 * j + 1]);
            }
        }

        // Print the coarse residual grid (add this part)
        std::cout << "Coarse Residual Grid (Cycle " << cycle + 1 << "):" << std::endl;
        printGrid();

        // Solve the coarse grid problem
        MultigridSolver coarseSolver(coarseSize);
        coarseSolver.setGrid(coarseResidual); // Initialize the coarse grid with the residual
        coarseSolver.solve();

        // Prolongation
        double** correction = new double* [N];
        for (int i = 0; i < N; ++i) {
            correction[i] = new double[N]();
        }
        for (int i = 0; i < coarseSize; ++i) {
            for (int j = 0; j < coarseSize; ++j) {
                correction[2 * i][2 * j] = coarseSolver.getSolution(i, j);
                correction[2 * i + 1][2 * j] = coarseSolver.getSolution(i, j);
                correction[2 * i][2 * j + 1] = coarseSolver.getSolution(i, j);
                correction[2 * i + 1][2 * j + 1] = coarseSolver.getSolution(i, j);
            }
        }

        // Print the correction grid (add this part)
        std::cout << "Correction Grid (Cycle " << cycle + 1 << "):" << std::endl;
        printGrid();

        // Correction
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                u[i][j] += correction[i][j];
            }
        }

        // Post-smoothing
        smooth();

        // Print intermediate grid (add this part)
        std::cout << "Intermediate Grid (Cycle " << cycle + 1 << "):" << std::endl;
        printGrid();

        // Deallocate memory for the residual grid
        for (int i = 0; i < N; ++i) {
            delete[] residual[i];
        }
        delete[] residual;

        // Deallocate memory for the correction grid
        for (int i = 0; i < N; ++i) {
            delete[] correction[i];
        }
        delete[] correction;

        // Deallocate memory for the coarse residual grid
        for (int i = 0; i < coarseSize; ++i) {
            delete[] coarseResidual[i];
        }
        delete[] coarseResidual;
    }
}

