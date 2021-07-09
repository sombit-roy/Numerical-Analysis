# Numerical Analysis

Codes for solving definite integrals, ODEs and PDEs using numerical analysis techniques and displaying the physical systems described by them.

## Definite Integrals

In the Simpson method, a quadratic polynomial is interpolated from sets of three consecutive points (of spacing h) and the areas under them are added up. The formula to approximate definite integrals by this method is

<img src="https://render.githubusercontent.com/render/math?math=\int_{x_1}^{x_N} f(x)dx = \cfrac{h}{3}(f(x_1) %2B 4f(x_2) %2B 2f(x_3) %2B 4f(x_4) %2B 2f(x_5) %2B \cdots %2B 2f(x_{N-2}) %2B 4f(x_{N-1}) %2B f(x_N))">

We can get the diffraction pattern from the parametric curve known as Cornu's spiral

<img src="https://render.githubusercontent.com/render/math?math=C(t) = \int_0^t \cos\left(\frac{\pi u^2}{2}\right) du"><br>
<img src="https://render.githubusercontent.com/render/math?math=S(t) = \int_0^t \sin\left(\frac{\pi u^2}{2}\right) du">

with <img src="https://render.githubusercontent.com/render/math?math=t \in [-\infty,\infty]">.

The intensity due to diffraction at a sharp edge is plotted as

<img src="https://render.githubusercontent.com/render/math?math=I(t) = {\left(\frac{1}{2}-C(t)\right)}^2 %2B {\left(\frac{1}{2}-S(t)\right)}^2">

![](figures/int_diffraction.png)

Another way of integration is the Monte Carlo method. The surface to be integrated is enclosed in a box. A large number of points are scattered at random. The proportion of points inside the surface to the total is equal to the proprtion of the area of the surface to that of the box.

The dipole potential of a disk lying in the x-y plane, with one half of positive charge and another half of negative charge, is

<img src="https://render.githubusercontent.com/render/math?math=V(\textbf{r}) = \cfrac{1}{4\pi\epsilon_0 r^2} \int r' \cos\theta \rho(\textbf{r}) d\tau'">

where

<img src="https://render.githubusercontent.com/render/math?math=d\tau'"> = infinitesimal volume of charged surface <br>
<img src="https://render.githubusercontent.com/render/math?math=r"> = distance of infinitesimal volume from origin <br>
<img src="https://render.githubusercontent.com/render/math?math=r'"> = distance of point of measurement from origin <br>
<img src="https://render.githubusercontent.com/render/math?math=\theta"> = angle between above two <br>
<img src="https://render.githubusercontent.com/render/math?math=\rho"> = charge distribution function

Using MC integration, we get the following plot for potential on the x-y plane

![](figures/int_dipole.png)

## Ordinary Differential Equations

The Runge Kutta method of order 4 can be used to solve an initial value ODE of the form <img src="https://render.githubusercontent.com/render/math?math=\cfrac{dy}{dx} = f(x,y)"> by constructing

<img src="https://render.githubusercontent.com/render/math?math=k_1 = f(x_n,y_n)"><br>
<img src="https://render.githubusercontent.com/render/math?math=k_1 = f(x_n %2B \frac{h}{2},y_n %2B \frac{hk_1}{2})"><br>
<img src="https://render.githubusercontent.com/render/math?math=k_1 = f(x_n %2B \frac{h}{2},y_n %2B \frac{hk_2}{2})"><br>
<img src="https://render.githubusercontent.com/render/math?math=k_1 = f(x_n %2B h,y_n %2B hk_3)">

We get <img src="https://render.githubusercontent.com/render/math?math=y_{n %2B 1} = y_n %2B \cfrac{1}{6}(k_1 %2B 2k_2 %2B 2k_3 %2B k_4) %2B \mathcal{O}(h^5)"> and solve this iteratively.

The two body problem is a coupled second order system. Its governing differential equation is Newton's gravitational law

<img src="https://render.githubusercontent.com/render/math?math=\cfrac{d^2 \textbf{r}_1}{d t^2} = \cfrac{-Gm_2 \hat{\textbf{r}}_{12}}{|\textbf{r}_{12}|^2} = \cfrac{-Gm_2 (\textbf{r}_1 - \textbf{r}_2)}{|\textbf{r}_1-\textbf{r}_2|^3}"><br>
<img src="https://render.githubusercontent.com/render/math?math=\cfrac{d^2 \textbf{r}_2}{d t^2} = \cfrac{-Gm_1 \hat{\textbf{r}}_{21}}{|\textbf{r}_{21}|^2} = \cfrac{-Gm_1 (\textbf{r}_2 - \textbf{r}_1)}{|\textbf{r}_2-\textbf{r}_1|^3}"><br>

We can decouple the equations by considering <img src="https://render.githubusercontent.com/render/math?math=y_1 = \textbf{r}_1, y_2 = \textbf{r}_2, y_3 = \dot{\textbf{r}_1}, y_4 = \dot{\textbf{r}_2}">. Thus, we can see that

<img src="https://render.githubusercontent.com/render/math?math=\dot{y_1} = y_3"><br>
<img src="https://render.githubusercontent.com/render/math?math=\dot{y_2} = y_4"><br>
<img src="https://render.githubusercontent.com/render/math?math=\dot{y_3} = \cfrac{-Gm_2 (y_1 - y_2)}{\left(\sqrt{(y_{1x} - y_{2x})^2 %2B (y_{1y} - y_{2y})^2 %2B (y_{1z} - y_{2z})^2}\right)^3}"><br>
<img src="https://render.githubusercontent.com/render/math?math=\dot{y_4} = \cfrac{-Gm_1 (y_2 - y_1)}{\left(\sqrt{(y_{1x} - y_{2x})^2 %2B (y_{1y} - y_{2y})^2 %2B (y_{1z} - y_{2z})^2}\right)^3}">

We solve them using the RK method with <img src="https://render.githubusercontent.com/render/math?math=k_{11} = f_1(x_n,y_{1n},y_{2n}), k_21 = f_2(x_n,y_{1n},y_{2n})">, and so on for <img src="https://render.githubusercontent.com/render/math?math=k_{12}, k_{22}">, etc. The function returns updated position and velocity vectors.

As an example, initially we take the positions as <img src="https://render.githubusercontent.com/render/math?math=\textbf{r}_1 = 10\hat{x}, \textbf{r}_2 = -10\hat{x}"> and their corresponding velocities as <img src="https://render.githubusercontent.com/render/math?math=\textbf{v}_1 = 0.1\hat{y} %2B 0.1\hat{z}, \textbf{v}_2 = -0.1\hat{y}">. (Note: We take G = 1)

Below plotted are the trajectories of the particles for <img src="https://render.githubusercontent.com/render/math?math=m_1 = 1.1, m_2 = 0.8">

![](figures/ode_gravity_1.gif)

We again plot for <img src="https://render.githubusercontent.com/render/math?math=m_1 = 1.1, m_2 = 0.008">. As expected, if one of the masses is much heavier than the other, the smaller particle revolves around the bigger one.

![](figures/ode_gravity_2.gif)

For the ODE <img src="https://render.githubusercontent.com/render/math?math=y''=f(x,y,y')">, given <img src="https://render.githubusercontent.com/render/math?math=a<x<b">, the boundary conditions are defined as <img src="https://render.githubusercontent.com/render/math?math=y(a)=y_a"> and <img src="https://render.githubusercontent.com/render/math?math=y(b)=y_b">. The shooting method is used to convert the boundary value problem to an initial value problem with initial conditions <img src="https://render.githubusercontent.com/render/math?math=y(a)=y_a"> and <img src="https://render.githubusercontent.com/render/math?math=y'(a)=\alpha^{(k)}">, where <img src="https://render.githubusercontent.com/render/math?math=\alpha^{(k)}"> is such that it satisfies <img src="https://render.githubusercontent.com/render/math?math=y(b)=y_b">.

<img src="https://render.githubusercontent.com/render/math?math=\alpha^{(0)}"> is initially chosen as the slope between points <img src="https://render.githubusercontent.com/render/math?math=(a,y_a)"> and <img src="https://render.githubusercontent.com/render/math?math=(b,y_b)">. This trial <img src="https://render.githubusercontent.com/render/math?math=\alpha"> is passed as the argument <img src="https://render.githubusercontent.com/render/math?math=y_b"> to the Runge-Kutta function, which returns a vector <img src="https://render.githubusercontent.com/render/math?math=y"> whose last element is an updated value of <img src="https://render.githubusercontent.com/render/math?math=y_b">.

Now, we must solve for <img src="https://render.githubusercontent.com/render/math?math=\alpha^{(k)}"> iteratively until <img src="https://render.githubusercontent.com/render/math?math=y_b^{(k)}"> is within the margin of error of the actual <img src="https://render.githubusercontent.com/render/math?math=y_b">. The recurrence relation is given by comparing slopes

<img src="https://render.githubusercontent.com/render/math?math=\alpha^{(k)} = \alpha^{(k-2)} %2B (y_b - y_b^{(k-2)})\cfrac{\alpha^{(k-1)} - \alpha^{(k-2)}}{y_b^{(k-1)} - y_b^{(k-2)}}">

The hydrogen wavefunction is <img src="https://render.githubusercontent.com/render/math?math=\Psi(r,\theta,\phi) = R(r)\Theta(\theta)\Phi(\phi)">. The azimuthal shape of the orbital with quantum numbers l and m is given by <img src="https://render.githubusercontent.com/render/math?math=\Theta(\theta) = P_{lm}(\cos\theta)">, which is the solution of the associated Legendre equation 

<img src="https://render.githubusercontent.com/render/math?math=(1-x^2)\cfrac{d^2}{d x^2} P_{lm}(x) - 2x\cfrac{d}{dx} P_{lm}(x) %2B \left(l(l %2B 1) - \cfrac{m^2}{1-x^2}\right) P_{lm}(x) = 0">

with boundary conditions <img src="https://render.githubusercontent.com/render/math?math=y(-1)=(-1)^n"> and <img src="https://render.githubusercontent.com/render/math?math=y(1)=1">.

We map <img src="https://render.githubusercontent.com/render/math?math=x \mapsto \cos\theta"> plot the square of the associated Legendre polynomial to get the electron density distribution.

This is the orbital shape of <img src="https://render.githubusercontent.com/render/math?math=p_z">, i.e. l=1, m=0.

![](figures/ode_orbital_1.png)

And this is the orbital shape of <img src="https://render.githubusercontent.com/render/math?math=d_{xz}">, i.e. l=2, m=1.

![](figures/ode_orbital_2.png)

## Partial Differential equations

To solve PDEs we use the finite difference method, in which spatial and temporal discretization is done.

We use this method to solve the heat equation for a spherical body. A metal ball at room temperature (25 degrees Celcius) is immersed in boiling water (100 degrees Celsius) till the temperature at the center is 40 degrees Celsius. The equation governing heat flow is given by

<img src="https://render.githubusercontent.com/render/math?math=\cfrac{1}{r^2}\cfrac{\partial}{\partial r}\left(r^2 \cfrac{\partial{T}}{\partial r}\right) = \cfrac{\rho c}{K} \cfrac{\partial T}{\partial t} \implies \cfrac{\partial T}{\partial t} = a\left(\cfrac{\partial ^2 T}{\partial r^2} %2B \cfrac{2}{r}\cfrac{\partial T}{\partial r}\right)"> where <img src="https://render.githubusercontent.com/render/math?math=a = \cfrac{K}{\rho c}">

Initially <img src="https://render.githubusercontent.com/render/math?math=T_{\forall r,t=0} = T_{\text{room}}"> and at the boundary <img src="https://render.githubusercontent.com/render/math?math=T_{r=R,\forall t} = T_{\text{water}}">

Now, after applying the finite difference method, we get

<img src="https://render.githubusercontent.com/render/math?math=\cfrac{T_{r,t %2B 1} - T_{r,t}}{ak} = \cfrac{2}{r}\left(\cfrac{T_{r %2B 1,t} - T_{r-1,t}}{2h}\right) %2B \cfrac{T_{r %2B 1,t} - 2T_{r,t} %2B T_{r-1,t}}{h^2} \implies T_{r,t %2B 1} = (1-2\alpha)T_{r,t} %2B \cfrac{\alpha}{r}\left((r-1)T_{r-1,t} %2B (r %2B 1)T_{r %2B 1,t} \right)"> where <img src="https://render.githubusercontent.com/render/math?math=\alpha = \cfrac{ak}{h^2}">

![](figures/pde_heat.gif)

Using the finite difference method, we can also solve the wave equation for a string with a mass bead on it (located at three-fourths of its length from the left end) and show the reflection and tranmission of waves caused by it.

<img src="https://render.githubusercontent.com/render/math?math=\cfrac{\partial ^2 y}{\partial x^2} = \cfrac{1}{c^2}\cfrac{\partial ^2 y}{\partial c^2}">

The initial wave is of Gaussian nature given by <img src="https://render.githubusercontent.com/render/math?math=y_{x\in(-1,1),t=0} = e^{\frac{-1}{1-x^2}}">. It is given a right moving velocity as <img src="https://render.githubusercontent.com/render/math?math=\cfrac{1}{c}\cfrac{\partial y}{\partial x} = \cfrac{2xy}{(1-x^2)^2}">.

We discretize the equation, with an extra condition relating to the mass bead 

<img src="https://render.githubusercontent.com/render/math?math=m_{\text{bead}} \left(\cfrac{y_{x,t %2B 1} - 2y_{x,t} %2B y_{x,t-1}}{k^2}\right) = -T\sin\theta_1 %2B T\sin\theta_2 \approx -T\tan\theta_1 %2B T\tan\theta_2 = -T\left(\cfrac{y_{x,t} - y_{x-1,t}}{h}\right)%2B T\left(\cfrac{y_{x %2B 1,t} - y_{x,t}}{h}\right)">

![](figures/pde_string.gif)