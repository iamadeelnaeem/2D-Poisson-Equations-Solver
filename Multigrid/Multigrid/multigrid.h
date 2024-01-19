// multigrid.h

#ifndef MULTIGRID_H
#define MULTIGRID_H

class MultigridSolver {
public:
    MultigridSolver(int size);
    ~MultigridSolver();

    void initializeRHS();
    void solve();
    void setGrid(double** initialGuess);
    double getSolution(int i, int j);
    void setBoundaryValues();
    void initializeGrid();
    void printGrid() const;

private:
    int N; // Size of the grid
    double** u; // Solution grid
    double** rhs; // Right-hand side
    // Add any other necessary member variables

    void smooth();
    void restrict();
    void prolongate();
    
};

#endif // MULTIGRID_H
