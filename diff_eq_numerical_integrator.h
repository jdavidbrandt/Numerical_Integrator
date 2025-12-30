#ifndef DIFF_EQ_NUMERICAL_INTEGRATOR_H_INCLUDED
#define DIFF_EQ_NUMERICAL_INTEGRATOR_H_INCLUDED

#include <cmath>
#include <vector>
#include <iostream>
#include <functional>
#include <cstdarg>



//takes arrays, with only 2 numbers, the first number in the array being the independent, the second being the dependent
class InitialCondition{
    public:
        
        //The initial condition y(pi) = 5 in this array would read, initial_condition[0] == pi and initial_condition[1] == 5, pairing the point (pi, 5) 
        double pairs[2];

        //constructor
        InitialCondition(double passed_independent_0, double passed_dependent_0){
            pairs[0] = passed_independent_0;
            pairs[1] = passed_dependent_0;
        }

};

//base class
class DiffEqn {
    public:
        //unable to make equation part of base class, because as order goes up, more doubles are required as input
        std::vector<double> predicted_rk4_independent_vals;
        std::vector<double> predicted_rk4_dependent_vals;
};

//of form dy/dt + P(t)y = Q(t). Define Diff eqn by solving for dy/dt, and returning what is on the right hand side of the eqn in main. dy/dt = return value
class FirstOrderODE: public DiffEqn{
    public:
        //define initial conditions, initiate a function (empty, it is defined in main), and initiate vectors for independent and dependent predicted values
        std::function<double(double,double)> equation;
        double independent_initial;
        double dependent_initial;

        //takes in what dy/dt is equal to (passed_equation), and an object of type InitialCondition
        FirstOrderODE(std::function<double(double,double)> passed_equation, InitialCondition initial_condition){
            equation = passed_equation;
            independent_initial = initial_condition.pairs[0];
            dependent_initial = initial_condition.pairs[1];

            //sets the first step
            predicted_rk4_independent_vals.push_back(initial_condition.pairs[0]);
            predicted_rk4_dependent_vals.push_back(initial_condition.pairs[1]);

            
        }

    //function to solve the the First Order Linear Diff Eqn, given a step size, and independent axis range
        void FirstOrderRK4Solve(double step_size, double independent_final){

            //Runge-Kutta constants initialization, and setting the first step to the initial conditions
            double k1, k2, k3, k4;
            double independent_step = independent_initial;
            double dependent_step = dependent_initial;
            

            //Define number of steps, and start Estimating slopes at every step, give them a weighted average, and fill vectors so that one can plot predicted values
            int n = (independent_final - independent_initial)/step_size;
            for (int i=0; i<n; ++i){

                //if independent_step == 0, and we are dividing by the current independent step, go to the next step, and don't calculate k vals
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

                    //adds each predicted point to two vectors, such that one can plot.
                    predicted_rk4_independent_vals.push_back(independent_step);
                    predicted_rk4_dependent_vals.push_back(dependent_step);
                }   
            }
        }
};

class SecondOrderODE: public DiffEqn{
    public:
        //populate constructor, then split into two first order ODE's and solve them using the FirstOrderRK4Solve method
        std::function<double(double,double,double)> equation;
        std::function<double(double,double,double)> v_dot; 
        double independent_initial, dependent_initial, independent_prime_initial, dependent_prime_initial;
        std::vector<double> predicted_prime_rk4_independent_vals, predicted_prime_rk4_dependent_vals;

        //takes in what d2y/dt2 is equal to, i.e. d2y/dt2 = f(t, y, y'), and a vector populated with types InitialCondtion (there should be two Initial conditions in this particular vector)
        SecondOrderODE(std::function<double(double,double,double)> passed_equation, std::vector<InitialCondition> vInitialConditions){
            //sets equations, and initial conditions. The equations are immediately broken into two first order differential equations
            //where v = f1(t, y, y_dot) and v_dot = f2(t, y, y_dot) 

            //split second order diff eqn into two first order by defining two fundtions, one for v_dot (so copy the passed function), and one for v, i.e, x_dot
            equation, v_dot = passed_equation;

            independent_initial = vInitialConditions.at(0).pairs[0];
            dependent_initial = vInitialConditions.at(0).pairs[1];
            independent_prime_initial = vInitialConditions.at(1).pairs[0];
            dependent_prime_initial = vInitialConditions.at(1).pairs[1];

            //sets the first steps
            predicted_rk4_independent_vals.push_back(independent_initial);
            predicted_rk4_dependent_vals.push_back(dependent_initial);
            predicted_prime_rk4_independent_vals.push_back(independent_prime_initial);
            predicted_prime_rk4_dependent_vals.push_back(dependent_prime_initial);

            
            
        }

        //takes in step_size of type double, and independent_final of type double
        void SecondOrderRK4Solve(double step_size, double independent_final){
            //initialize Runge-Kutta constants k being for dependent, and l being for dependent_prime, along with the substitution equation, v
            std::function<double(double,double,double)> v = [](double t, double y, double y_dot){
                return y_dot;
            };
            double k1, k2, k3, k4, l1, l2, l3, l4 = 0;
            double independent_step = independent_initial;
            double dependent_step = dependent_initial;
            double independent_prime_step = independent_prime_initial;
            double dependent_prime_step = dependent_prime_initial;

            int n = (independent_final - independent_initial)/step_size;

            for (int i = 0; i < n; ++i){
                if(std::isnan(v_dot(independent_step, dependent_step, dependent_prime_step)) || std::isinf(v_dot(independent_step, dependent_step, dependent_prime_step))){
                    independent_step+=step_size;
                }
                else{
                    k1 = step_size*v(independent_step, dependent_step, dependent_prime_step);
                    l1 = step_size*v_dot(independent_step, dependent_step, dependent_prime_step);
                    k2 = step_size*v(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k1, dependent_prime_step + 0.5*step_size*l1);
                    l2 = step_size*v_dot(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k1, dependent_prime_step + 0.5*step_size*l1);
                    k3 = step_size*v(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k2, dependent_prime_step + 0.5*step_size*l2);
                    l3 = step_size*v_dot(independent_step + 0.5*step_size, dependent_step + 0.5*step_size*k2, dependent_prime_step + 0.5*step_size*l2);
                    k4 = step_size*v(independent_step + step_size, dependent_step + step_size*k3, dependent_prime_step + step_size*l3);
                    l4 = step_size*v_dot(independent_step + step_size, dependent_step + step_size*k3, dependent_prime_step + step_size*l3);

                    independent_step += step_size;
                    independent_prime_step += step_size;
                    dependent_step += (k1 + 2*k2 + 2*k3 + k4)/6.0;
                    dependent_prime_step += (l1 + 2*l2 + 2*l3 + l4)/6.0;

                    predicted_rk4_independent_vals.push_back(independent_step);
                    predicted_rk4_dependent_vals.push_back(dependent_step);
                    predicted_prime_rk4_independent_vals.push_back(independent_prime_step);
                    predicted_prime_rk4_dependent_vals.push_back(dependent_prime_step);
                }
            }
        }        
        

};

class HigherOrderODE: public DiffEqn{
    public:
        //takes the order of the differential equation, the passed equation (with coefficent of 1 for highest order), and a vector of initial conditions sorted from 0th derivative to nth order derivative
        int order;
        std::vector<InitialCondition> vInitialConditions;

        HigherOrderODE(int passed_order, std::function<double(double,double)> passed_equation, std::vector<InitialCondition> passed_InitialConditions){
            order = passed_order;
            //equation = passed_equation;
            vInitialConditions = passed_InitialConditions;

        }
    
};























#endif