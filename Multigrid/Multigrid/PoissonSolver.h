// PoissonSolver.h
#pragma once

#include <vector>

class PoissonSolver {
public:
    PoissonSolver(int gridSize);

    void initializeGrid();
    void setBoundaryValues();
    void solve();
    void printSolution();

private:
    int gridSize;
    std::vector<std::vector<double>> grid;
    std::vector<std::vector<double>> newGrid;  // Temporary grid for Jacobi method
    int maxIterations;
    double tolerance;

    void jacobiSmoothing();
    void restrict();
    void interpolate();
};
