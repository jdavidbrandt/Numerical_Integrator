#ifndef DIFF_EQ_NUMERICAL_INTEGRATOR_H_INCLUDED
#define DIFF_EQ_NUMERICAL_INTEGRATOR_H_INCLUDED

#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <cstdarg>
#include <tuple>

/////////////////////////////////////////////////////////////////////////////////

//DEFINE OBJECT OF TYPE INITIAL CONDITION, SO PASSING THESE STATEMENTS IS EASIER

/////////////////////////////////////////////////////////////////////////////////

//takes arrays, with only 2 numbers, the first number in the array being the independent, the second being the dependent
class InitialCondition{
public:
    int order;
    double independent_0;
    double dependent_0;

    //The initial condition y(pi) = 5 in this array would read, initial_condition[0] == pi and initial_condition[1] == 5, pairing the point (pi, 5) 
    double initial_condition[2] = {0,0};

    //constructor
    InitialCondition(double passed_independent_0, double passed_dependent_0){
        initial_condition[0] = passed_independent_0;
        initial_condition[1] = passed_dependent_0;
    }

};

//of form dy/dt + P(t)y = Q(t). Define Diff eqn by solving for dy/dt, and returning what is on the right hand side of the eqn in main. dy/dt = return value
class FirstOrderODE{
public:
    //define initial conditions, initiate a function (empty, it is defined in main), and set vectors for independent and dependent predicted values
    double independent_initial;
    double dependent_initial;
    std::function<double(double, double)> equation;
    std::vector<double> predicted_rk4_independent_vals = {independent_initial};
    std::vector<double> predicted_rk4_dependent_vals = {dependent_initial};

    //constructor
    FirstOrderODE(std::function<double(double,double)> passed_equation, double passed_independent_0, double passed_dependent_0){
        equation = passed_equation;
        independent_initial = passed_independent_0;
        dependent_initial = passed_dependent_0;
    }

    //function to solve the the First Order Linear Diff Eqn, given a step size, and independent axis range
    void FirstOrderRK4Solve(double step_size, double independent_final){

        //Runge-Kutta constants initialization, and setting the first step to the initial conditions
        double k1, k2, k3, k4;
        double dependent_step = dependent_initial;
        double independent_step = independent_initial;

        //Define number of steps, and start Estimating slopes at every step, give them a weighted average, and fill vectors so that one can plot predicted values
        int n = (independent_final - independent_initial)/step_size;
        for (int i=0; i<n; ++i){

            //if independent_step == 0, go to next step, and don't calculate k vals
            if(std::isnan(equation(independent_step, dependent_step)) || std::isinf(equation(independent_step, dependent_step))){
                independent_step+=step_size;
            }
            else{
                k1 = step_size * equation(independent_step, dependent_step);
                k2 = step_size * equation(independent_step + 0.5*step_size, dependent_step + 0.5*k1);
                k3 = step_size * equation(independent_step + 0.5*step_size, dependent_step + 0.5*k2);
                k4 = step_size * equation(independent_step + step_size, dependent_step + k3);
                independent_step += step_size;
                dependent_step += (k1 + 2*k2 + 2*k3 + k4)/6.0;

                predicted_rk4_independent_vals.push_back(independent_step);
                predicted_rk4_dependent_vals.push_back(dependent_step);
            }
            
        }

    }
};























#endif