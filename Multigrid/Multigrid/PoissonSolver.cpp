// PoissonSolver.cpp
#include "PoissonSolver.h"
#include <iostream>

PoissonSolver::PoissonSolver(int gridSize) : gridSize(gridSize) {
    maxIterations = 100;
    tolerance = 1e-6;
}

void PoissonSolver::initializeGrid() {
    grid.resize(gridSize, std::vector<double>(gridSize, 0.0));
    newGrid.resize(gridSize, std::vector<double>(gridSize, 0.0));  // Initialize newGrid
}

void PoissonSolver::setBoundaryValues() {
    // Ask the user to input grid points on the boundary
    std::cout << "Enter a singal value for all the grid points on the boundary:" << std::endl;
    int ip;
    std::cin >> ip;
    for (int i = 0; i < gridSize; ++i) {

        grid[i][0] = ip;

        grid[i][gridSize - 1] = ip;
    }
    for (int j = 1; j < gridSize - 1; ++j) {
        grid[0][j] = ip;
        grid[gridSize - 1][j] = ip;
    }
}

void PoissonSolver::solve() {
    double maxChange;  // Maximum change in the solution
    double currentChange;  // Change in the current iteration

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        jacobiSmoothing();

        // Compute the maximum change in the solution
        maxChange = 0.0;
        for (int i = 1; i < gridSize - 1; ++i) {
            for (int j = 1; j < gridSize - 1; ++j) {
                currentChange = std::abs(grid[i][j] - newGrid[i][j]);
                if (currentChange > maxChange) {
                    maxChange = currentChange;
                }
            }
        }

        // Check for convergence
        if (maxChange < tolerance) {
            std::cout << "Converged after " << iteration + 1 << " iterations.\n";
            break;
        }
    }
}


void PoissonSolver::printSolution() {
    std::cout << "Solved Grid:\n";
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            std::cout << grid[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void PoissonSolver::jacobiSmoothing() {
    newGrid = grid;  // Initialize newGrid with the current values of grid

    for (int i = 1; i < gridSize - 1; ++i) {
        for (int j = 1; j < gridSize - 1; ++j) {
            newGrid[i][j] = 0.25 * (grid[i - 1][j] + grid[i + 1][j] + grid[i][j - 1] + grid[i][j + 1]);
        }
    }

    grid = newGrid;  // Update grid with the values from newGrid
}

void PoissonSolver::restrict() {
    int coarseSize = (gridSize - 1) / 2 + 1;  // Size of the coarser grid

    std::vector<std::vector<double>> coarseGrid(coarseSize, std::vector<double>(coarseSize, 0.0));

    // Perform restriction (average values)
    for (int i = 1; i < coarseSize - 1; ++i) {
        for (int j = 1; j < coarseSize - 1; ++j) {
            coarseGrid[i][j] = 0.25 * (grid[2 * i][2 * j] + grid[2 * i + 1][2 * j] +
                grid[2 * i][2 * j + 1] + grid[2 * i + 1][2 * j + 1]);
        }
    }

    // Update the grid with the coarser grid
    grid = coarseGrid;
    gridSize = coarseSize;
}


void PoissonSolver::interpolate() {
    int fineSize = 2 * gridSize - 1;  // Size of the finer grid

    std::vector<std::vector<double>> fineGrid(fineSize, std::vector<double>(fineSize, 0.0));

    // Perform interpolation using linear interpolation
    for (int i = 0; i < fineSize; ++i) {
        for (int j = 0; j < fineSize; ++j) {
            if (i % 2 == 0 && j % 2 == 0) {
                fineGrid[i][j] = grid[i / 2][j / 2];
            }
            else if (i % 2 == 0) {
                fineGrid[i][j] = 0.5 * (grid[i / 2][j / 2 - 1] + grid[i / 2][j / 2]);
            }
            else if (j % 2 == 0) {
                fineGrid[i][j] = 0.5 * (grid[i / 2 - 1][j / 2] + grid[i / 2][j / 2]);
            }
            else {
                fineGrid[i][j] = 0.25 * (grid[i / 2 - 1][j / 2 - 1] + grid[i / 2 - 1][j / 2 + 1] +
                    grid[i / 2 + 1][j / 2 - 1] + grid[i / 2 + 1][j / 2 + 1]);
            }
        }
    }

    // Update the grid with the finer grid
    grid = fineGrid;
    gridSize = fineSize;
}

