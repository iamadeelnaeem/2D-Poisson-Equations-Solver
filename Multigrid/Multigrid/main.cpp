// main.cpp
#include "PoissonSolver.h"
#include <iostream>

int main() {
    int gridSize;
    std::cout << "Enter grid size: ";
    std::cin >> gridSize;

    PoissonSolver solver(gridSize);
    solver.initializeGrid();
    solver.setBoundaryValues();
    solver.solve();
    solver.printSolution();

    return 0;
}
