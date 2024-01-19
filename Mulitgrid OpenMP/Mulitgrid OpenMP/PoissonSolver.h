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
    void restrict();
    void interpolate();

private:
    int gridSize;
    std::vector<std::vector<double>> grid;
    std::vector<std::vector<double>> newGrid;  // Temporary grid for Jacobi method
    int maxIterations;
    double tolerance;

    void jacobiSmoothing();
    
};