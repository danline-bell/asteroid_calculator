function [x,f,ea,iter]=newtmult(func,es,maxit,varargin) 
% newtmult: Newton-Raphson root zeroes nonlinear systems 
%  [x,f,ea,iter] = newtmult(func,x0,es,maxit,p1,p2,...): 
%    uses the Newton-Raphson method to find the roots of 
%    a system of nonlinear equations 
% input: 
%   func = name of function that returns f and J 
%   x0 = initial guess 
%   es = desired percent relative error (default = 0.0001%) 
%   maxit = maximum allowable iterations (default = 50) 
%   p1,p2,... = additional parameters used by function 
% output: 
%   x = vector of roots 
%   f = vector of functions evaluated at roots 
%   ea = approximate percent relative error (%) 
%   iter = number of iterations
% global variables
x0 = [1.5,3.5];
g = 9.81;
r_earth = 6371;
r_platform = 500;
steel_strength = 975; % N/mm2

if nargin<1,error('at least 2 input arguments required'),end 
if nargin<3|isempty(es),es = 0.0001;end 
if nargin<4|isempty(maxit),maxit = 50;end 
iter = 0; 
x = x0; 
while (1)  
    [J,f] = func(x,varargin{:});  
    dx = J\f;  
    x = x - dx;  
    iter = iter + 1;  
    ea=100*max(abs(dx./x));  
    if iter>= maxit|ea<= es, break, end 
end

% asteroid plot
x_aster = 0:100:50000;
y_aster = (32000*x_aster-(16/25)*x_aster.^2).^(1/2);
figure(1);
plot(x_aster,y_aster, 'b');
xlim([-5000 55000])
grid on
axis equal
hold on 
plot(x_aster,-y_aster,'b');

% earth circumference
angle = 0:0.01:2*pi;
xp = r_earth*cos(angle);
yp = r_earth*sin(angle);
plot(1e4+xp,yp,'g') % plot at center of earth

% platform location
x_platform = 13.186e3;
y_platform = 5.517e3;
angle = 0:0.01:2*pi;
xp = r_platform*cos(angle);
yp = r_platform*sin(angle);
plot(x_platform+xp,y_platform+yp,'k')

% rocket path
x_rocket = 13186:10:21428;
y_rocket = 1.7321*x_rocket-17321;
plot(x_rocket,y_rocket,'--c')

% strike location
x_strike = 21428;
y_strike = 19795;
plot(x_strike, y_strike, 'xr');

% rocket distance
d_rocket = ((x_strike-x_platform)^2+(y_strike-y_platform)^2)^(1/2)

% weight calculations (kg-km-s)
m_rocket = 1000; % mass of rocket
w_rocket = 1000 * g; % weight of rocket without fuel
ro_fuel = 750; % density of fuel kg/m3
fuel_rate = 0.25; % fuel burn rate L/km
V_fuel = 0.001 * fuel_rate * d_rocket;
volume_fuel = V_fuel * 1000 % volume of fuel in m^3
% m = ro * V, so
m_fuel = ro_fuel * V_fuel
w_fuel = m_fuel * g;
% total weight 
w_total = w_rocket + w_fuel

% truss members axial forces
% AB AC BC BD CD CE DE Ax Ay Ey
M = zeros(10,10);
M(1,1) = -cos(pi/3);
M(1,2) = 1;
M(1,8) = 1;
M(2,1) = -sin(pi/3);
M(2,9) = 1;
M(3,1) = cos(pi/3);
M(3,3) = -cos(pi/3);
M(3,4) = -1;
M(4,1) = sin(pi/3);
M(4,3) = sin(pi/3);
M(5,2) = -1;
M(5,3) = cos(pi/3);
M(5,5) = -cos(pi/3);
M(5,6) = 1;
M(6,3) = -sin(pi/3);
M(6,5) = -sin(pi/3);
M(7,4) = 1;
M(7,5) = cos(pi/3);
M(7,7) = -cos(pi/3);
M(8,5) = sin(pi/3);
M(8,7) = sin(pi/3);
M(9,7) = cos(pi/3);
M(9,6) = -1;
M(10,7) = -sin(pi/3);
M(10,10) = 1;

b = [0;0;0;w_total/2;0;0;0;w_total/2;0;0];

% axial forces to solve for
axial_forces = M\b;

% maximum axial force
max_axial_force = max(axial_forces)

% cross section 
cross_section = (2 * max_axial_force)/steel_strength
end
function [J,f]=jfreact(x,varargin) 
del = 0.000001; 
df1dx1 = (u(x(1)+del*x(1),x(2))-u(x(1),x(2)))/(del*x(1)); 
df1dx2 = (u(x(1),x(2)+del*x(2))-u(x(1),x(2)))/(del*x(2)); 
df2dx1 = (v(x(1)+del*x(1),x(2))-v(x(1),x(2)))/(del*x(1)); 
df2dx2 = (v(x(1),x(2)+del*x(2))-v(x(1),x(2)))/(del*x(2)); 
J=[df1dx1 df1dx2;df2dx1 df2dx2]; 
f1=u(x(1),x(2)); 
f2=v(x(1),x(2)); 
f=[f1;f2];
end
function f = u(x,y) 
f = (32000*x-(16/25)*x^2)^(1/2) - y;
end
function f = v(x,y) 
f = 1.7321 * x - y - 17321;
end
