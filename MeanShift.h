#include <vector>

using namespace std;

class MeanShift {
public:
    MeanShift() { set_kernel(NULL); }
    MeanShift(double (*_kernel_func)(double,double)) { set_kernel(kernel_func); }
    vector<vector<double> > cluster(vector<vector<double> >, double);
private:
    double (*kernel_func)(double,double);
    void set_kernel(double (*_kernel_func)(double,double));
    vector<double> shift_point(const vector<double> &, const vector<vector<double> > &, double);
    //The following is to mimic the python version of the meanshift for the labeling
    vector<int> group_points(vector<vector<double> > &);
};
