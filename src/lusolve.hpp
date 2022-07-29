

#include <vector>

int ludcmp(std::vector<std::vector<double> > &a, int n, std::vector<int> &indx, int *d);
void lubksb(std::vector<std::vector<double> >  &a, int n, std::vector<int> &indx, std::vector<double> &b);
int lusolve(std::vector<std::vector<double> >  &a, int n, std::vector<double> &b);

