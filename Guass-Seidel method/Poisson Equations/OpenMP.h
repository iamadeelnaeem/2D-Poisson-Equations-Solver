#ifndef OpenMP_H
#define OpenMP_H

#include <vector>

class OpenMP {
public:
    OpenMP(int size);
    ~OpenMP();

    void solve(std::vector<double>& phi, const std::vector<double>& f, double tolerance, int maxIterations);
    void preSmooth(std::vector<double>& x, const std::vector<double>& b);
    void postSmooth(std::vector<double>& x, const std::vector<double>& b);
    // Other member functions...

private:
    void restriction(const std::vector<double>& r, std::vector<double>& r2h, int n);
    void interpolation(const std::vector<double>& E2h, std::vector<double>& E, int n);
    std::vector<double> cycle(const std::vector<double>& r, int lambda, int theta);
    
    // Other private member functions...

    int size;  // Size of the grid
};

#endif  // OpenMP_H
