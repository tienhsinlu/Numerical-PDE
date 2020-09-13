function [t, y] = leap_frog(f,tspan,y0,dt)
%
% [t,y] = rk4(f,tspan,y0,h)
%
% A simple integration routine to solve the
% initial value ODE   y' = f(t,y), y(a) = y0,
% using the clascical 4-stage Runge-Kutta method
% with a fixed step size h.
% tspan = [a b] is the integration interval.
% Note that y and f can be vector functions


m = length(y0);		     % problem size
t= tspan(1):dt:tspan(2);  % output abscissae
N = length(t)-1;	     % number of steps
y = zeros(m,N+1);

y0 = y0(:);
[~,y1] = rk4(f, [0, dt], y0, dt);% make sure y0 is a column vector
y1=y1(:, end);
y(:,1) = y0;% initialize
y(:,2) = y1;
  
% Integrate
for i=2:N
    y_temp = y0+2*dt*f(t(i), y1);
    y0=y1;
    y1=y_temp;
    y(:, i+1)=y1;
end