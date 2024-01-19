// main.cpp
#include "PoissonSolver.h"
#include <iostream>

int main() {
    int gridSize;

    // Get user input for grid size
    std::cout << "Enter grid size: ";
    std::cin >> gridSize;

    // Create PoissonSolver object
    PoissonSolver solver(gridSize);

    // Initialize grid and set boundary values
    solver.initializeGrid();
    solver.setBoundaryValues();

    // Solve the Poisson equation using multigrid with OpenMP
#pragma omp parallel
    {
        // Use OpenMP to parallelize the solver operations
#pragma omp sections
        {
#pragma omp section
            solver.solve();
        }
    }

    // Print the solved grid
    solver.printSolution();

    return 0;
}
