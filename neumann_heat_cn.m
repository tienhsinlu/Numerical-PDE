% [u,err,x,t] = neumann_heat_cn(t_0,t_f,M,N)
%
% this solves the heat equation u_t = u_xx with initial data u_0 =
% cos(x) with Neumann boundary conditions using finite-differences in
% space and Crank-Nicolson time-stepping.  t_0 is the initial time, t_f is the
% final time, N is the number of mesh-points, and M is the number of
% time steps.  err is the error. 

function [u,err,x,t] = heat_cn(t_0,t_f,M,N)

% define the mesh in space
dx = 2*pi/(N-1);
x = 0:dx:2*pi;
x = x';

% define the mesh in time
dt = (t_f-t_0)/M;
t = t_0:dt:t_f;

% define the ratio r
r = dt/dx^2 

for i=1:N
  u(i,1) = cos(x(i));
end 
err(:,1) = u(:,1) - exp(t_0-t(1))*cos(x);

% for internal points, have
%    u_new(j) = u_old(j) + r/2*(u_old(j+1)-2*u_old(j)+u_old(j-1))
% for the two end-points, have
%    u_new(1) = u_old(1) + r/2*(2 u_old(2)-2*u_old(1))
%    u_new(N) = u_old(N) + r/2*( -2*u_old(N)+2 u_old(N-1))     
% clearly the endpoints are redundant: u(1)= u(N) at all times.  I just
% kept them around for plotting convenience.  

% define the matrix A which has to be inverted at every time-step.
%   u_new(1) - u_old(1) = dt/h^2 (2 u_new(2)-2*u_new(1))
%   u_new(i) - u_old(i) = dt/h^2 (u_new(i+1)-2*u_new(i)+u_new(i-1)
%   u_new(N) - u_old(N) = dt/h^2 ( -2*u_new(N)+ 2 u_new(N-1))
A = zeros(N-1,N-1);
A(1,1) = 1+r;
A(1,2) = -r;
for i=2:N-1
  A(i,i-1) = -r/2;
  A(i,i) = 1+r;
  A(i,i+1) = -r/2;
end
A(N,N) = 1+r;
A(N,N-1) = -r;

for j=1:M
  RHS(1) = u(1,j) + r/2*(2*u(2,j)-2*u(1,j));
  for i=2:N-1
    RHS(i) = u(i,j) + r/2*(u(i+1,j)-2*u(i,j)+u(i-1,j));
  end
  RHS(N) = u(N,j) + r/2*(-2*u(N,j)+2*u(N-1,j));
  u(:,j+1) = tri_diag(A,RHS);  
  err(:,j+1) = u(:,j+1) - exp(t_0-t(j+1))*cos(x);
end

