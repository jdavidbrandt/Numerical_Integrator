#include <iostream>
#include <cmath>
#include <matplot/matplot.h>
#include <vector>
#include "diff_eq_numerical_integrator.h"
#include "matrix.h"




int main(){
    double pi = 3.14159265;
    
    //FIRST ORDER ODE TESTING
    //////////////////////////////////////////////////////////////////////////
    //define diff eqn here and pass it in to constructor, dy/dt = return value
    std::function<double(double,double)> my_eqn = [](double t, double y){
        return -y*std::tan(t) + 2*std::cos(t)*std::cos(t)*std::sin(t) - (1.0/std::cos(t));
        
    };

    //define conditions, and differential equation
    double plot_start = 0;
    double t_initial = pi/4.0;
    double y_initial = 3*std::sqrt(2);
    double t_final= 20;
    double step_size = 0.00001;

    InitialCondition initial_condition(t_initial, y_initial);
    FirstOrderODE my_diffeqn(my_eqn, initial_condition);
    my_diffeqn.FirstOrderRK4Solve(step_size, plot_start, t_final);

    //actual solution to this particular diff eqn, this is to compare!
    std::vector<double> t_real = matplot::linspace(plot_start, t_final);
    std::vector<double> y_real;
    double y_val, t_val, sub;
    for (double val : t_real){
        y_val = -((1.0/2.0)*std::cos(val)*std::cos(2*val)) - std::sin(val) + 7*std::cos(val);
        
        y_real.push_back(y_val);
    }
    //////////////////////////////////////////////////////////////////

    //Solve diff eqn, and plot, compare to the actual solution above
    auto fig_one = matplot::figure(true);
    matplot::plot(my_diffeqn.predicted_rk4_independent_vals, my_diffeqn.predicted_rk4_dependent_vals, "g", t_real, y_real, "b--o");
    //matplot::show(fig_one);
    matplot::save("first_order_ode.png");
    

    //SECOND ORDER ODE TESTING
    
    //Differential Equation that yields Bessel Solutions. This particular one is for when alpha == 0,     https://en.wikipedia.org/wiki/Bessel_function#Bessel_functions_of_the_first_kind:_J%CE%B1
    std::function<double(double,double,double)> my_secondorder_eqn = [](double t, double y, double y_prime){
        return (-y*(t*t)/(t*t)) - y_prime*(1.0/t);
    };

    //Define initial conditions and step size
    double so_t_initial = 1;
    double so_t_final=20;
    double bessel_initial_val = std::cyl_bessel_j(0, so_t_initial);
    double bessel_prime_initial_val = -std::cyl_bessel_j(1, so_t_initial);
    double so_step_size = 0.0001;
    InitialCondition posInitialCond(so_t_initial, bessel_initial_val);
    InitialCondition velInitialCond(so_t_initial, bessel_prime_initial_val);


    //define SecondOrderODE passing the Bessel Diff Eqn of the First Kind, Order 0, and solve, member predicted_rk4_independent_vals and predicted_rk4_dependent_vals will plot
    SecondOrderODE my_second_order_ODE(my_secondorder_eqn, posInitialCond, velInitialCond);
    my_second_order_ODE.SecondOrderRK4Solve(so_step_size, bessel_initial_val, so_t_final);

    //////////////////////////////////////////////////////////////////
    //actual solution to this particular diff eqn, this is to compare!
    std::vector<double> so_t_real = matplot::linspace(bessel_initial_val, so_t_final);
    std::vector<double> so_y_real;
    std::vector<double> so_y_prime_real;
    double so_y_val, so_y_prime_val, so_t_val;
    for (double so_val : so_t_real){
        so_y_val = std::cyl_bessel_j(0.0, so_val);
        so_y_prime_val = -std::cyl_bessel_j(1.0, so_val);
        so_y_real.push_back(so_y_val);
        so_y_prime_real.push_back(so_y_prime_val);
    }
    //////////////////////////////////////////////////////////////////

    auto fig_two = matplot::figure(true);

    //plotting RK4 predictions with actual analytic solutions
    
    matplot::plot(my_second_order_ODE.predicted_rk4_independent_vals, my_second_order_ODE.predicted_rk4_dependent_vals, "g",
        so_t_real, so_y_real, "b--o",
        my_second_order_ODE.predicted_prime_rk4_independent_vals, my_second_order_ODE.predicted_prime_rk4_dependent_vals, "r",
        so_t_real, so_y_prime_real, "b--o");
        

    //matplot::plot(t_real, y_real, "b--o", t_real, y_prime_real, "b--o");
    //matplot::show(fig_two);
    matplot::save("bessel_solution.png");
    return 0;
}
