#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <vector>
#include "diff_eq_numerical_integrator.h"
#include "matrix.h"




int main(){
    double pi = 3.14159265;
    //FIRST ORDER ODE TESTING

    /*
    //define diff eqn here and pass it in to constructor, dy/dt = return value
    std::function<double(double,double)> my_eqn = [](double t, double y){
        return (-3*t*y + 3*t)/(pow(t,2) +1);
        
    };

    //define conditions, and differential equation
    double t_initial = 0;
    double y_initial = 2;
    double t_final=6;
    double step_size = 0.00001;

    InitialCondition initial_condition(t_initial, y_initial);
    FirstOrderODE my_diffeqn(my_eqn, initial_condition);

    //////////////////////////////////////////////////////////////////
    //actual solution to this particular diff eqn, this is to compare!
    std::vector<double> t_real = matplot::linspace(t_initial, t_final);
    std::vector<double> y_real;
    double y_val, t_val, sub;
    for (double val : t_real){
        sub = pow(val,2) + 1;
        y_val = 1 + 1/(pow(sub, 3.0/2.0));
        y_real.push_back(y_val);
    }
    //////////////////////////////////////////////////////////////////

    //Solve diff eqn, and plot, compare to the actual solution above
    my_diffeqn.FirstOrderRK4Solve(step_size, t_final);
    matplot::plot(my_diffeqn.predicted_rk4_independent_vals, my_diffeqn.predicted_rk4_dependent_vals, "g", t_real, y_real, "b--o");
    matplot::show();
    */

    //SECOND ORDER ODE TESTING
    std::function<double(double,double,double)> my_secondorder_eqn = [](double t, double y, double y_prime){
        return 2*y_prime - 5*y + 2*exp(-t)*(sin(2*t) - 2*cos(2*t));
    };

    //define vector of initial conditions;
    std::vector<InitialCondition> vInitialConditions;
    double t_initial=0;
    double t_final=6;
    double step_size = 0.0001;

    //Initial Conditions y(0) == 0 and y'(0) == 2, place them into vector of type InitialCondition
    InitialCondition initialcond0(0,0);
    InitialCondition initialcond1(0, 2);
    vInitialConditions.push_back(initialcond0);
    vInitialConditions.push_back(initialcond1);

    //define SecondOrderODE passing the equation y'' = ________ and the vector of InitialConditions
    SecondOrderODE my_second_order_ODE(my_secondorder_eqn, vInitialConditions);

    //add an argument at the end of this, t_start, that selects where I want the first point to be plotted (not an initial condition, just the lowest t value I want to start the plot)
    //then create an additional loop in diff_eq_numerical_integrator.h that goes back and calculates all of those points for y and y prime.
    my_second_order_ODE.SecondOrderRK4Solve(step_size, t_final);

    //////////////////////////////////////////////////////////////////
    //actual solution to this particular diff eqn, this is to compare!
    std::vector<double> t_real = matplot::linspace(t_initial, t_final);
    std::vector<double> y_real;
    std::vector<double> y_prime_real;
    double y_val, y_prime_val, t_val;
    for (double val : t_real){
        y_val = cosh(val)*sin(2*val);
        y_prime_val = sinh(val)*sin(2*val) + 2*cosh(val)*cos(2*val);
        y_real.push_back(y_val);
        y_prime_real.push_back(y_prime_val);
    }
    //////////////////////////////////////////////////////////////////

    //plotting RK4 predictions with actual analytic solutions
    matplot::plot(my_second_order_ODE.predicted_rk4_independent_vals, my_second_order_ODE.predicted_rk4_dependent_vals, "g",
        t_real, y_real, "b--o",
        my_second_order_ODE.predicted_prime_rk4_independent_vals, my_second_order_ODE.predicted_prime_rk4_dependent_vals, "r",
        t_real, y_prime_real, "b--o");
    matplot::show();

    return 0;
}
