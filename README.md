# Numerical_Integrator
Numerical Integrator using Runge-Kutta Order 4.

The first example is a First Order ODE

$\cos{t} \frac{dy}{dt} + \sin{t} y = 2 \cos^3{t} \sin{t} - 1$ with $y(\frac{\pi}{4}) = 3 \sqrt{2}$

In main.cpp, one can see that the blue-dotted line is the actual solution, where as the green line is the Runge-Kutta 4 prediction

With Solution stated here, and pictured below

$- \frac{1}{2} \cos{t} \cos{2t} - \sin{t} + 7 \cos{t}$

<img width="560" height="420" alt="first_order_ode" src="https://github.com/user-attachments/assets/738d6243-9879-4879-8be2-9566e075b90b" />





The next example is a Second Order ODE, with Bessel solutions.

$t^{2} \frac{d^{2}}{dt^{2}} + t \frac{dy}{dt} + (t^{2} - \alpha^{2}) y = 0$ with $y(1) = J_{0}(1)$ and $y'(1) = - J_{1}(1)$ 

In this case, $\alpha = 0$ so this solution is a 0th Order of the First Kind Bessel Function.

The actual solution is in blue, and the predicted solution with the Runge-Kutta 4 algorithm is shown in green

In blue, again, is the actual solution. One can plot this in desmos with the following and compare to the graph below

$J(a,t) = \sum{m=0}^{100} \frac{(-1)^{m}}{m!(m+a)!} (\frac{x}{2})^{2m + a}$
$J(0,t) {t > 0}$
$J(1,t) {t > 0)$
