#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <vector>
#include "diff_eq_numerical_integrator.h"
#include "matrix.h"




int main(){
    //using namespace matplot;
    double pi = 3.14159265;
    
    //define diff eqn here and pass it in to constructor, dy/dt = return value
    std::function<double(double,double)> my_eqn = [](double t, double y){
        return (-3*t*y + 3*t)/(pow(t,2) +1);
    };

    //define conditions, and differential equation
    double t_initial = 0;
    double y_initial = 2;
    InitialCondition initial_condition(t_initial, y_initial);
    double t_final=6;
    double step_size = 0.00001;
    FirstOrderODE my_diffeqn(my_eqn, initial_condition);

    //actual solution to the diff eqn
    std::vector<double> t_real = matplot::linspace(t_initial, t_final);
    std::vector<double> y_real;
    double y_val, t_val, sub;
    for (double val : t_real){
        sub = pow(val,2) + 1;
        y_val = 1 + 1/(pow(sub, 3.0/2.0));
        y_real.push_back(y_val);
    }

    //Solve diff eqn, and plot dependent variable, compare to the actual solution above
    my_diffeqn.FirstOrderRK4Solve(step_size, t_final);
    matplot::plot(my_diffeqn.predicted_rk4_independent_vals, my_diffeqn.predicted_rk4_dependent_vals, "g", t_real, y_real, "b--o");
    matplot::show();

    return 0;
}
