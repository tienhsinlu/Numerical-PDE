function fdstencil(k,j)

% Compute stencil coefficients for finite difference approximation 
% of k'th order derivative on a uniform grid.  Print the stencil and
% dominant terms in the truncation error.
%
% j should be a vector of indices of grid points, values u(x0 + j*h)
% are used, where x0 is an arbitrary grid point and h the mesh spacing.
% This routine returns a vector c of length n=length(j) and the
% k'th derivative is approximated by
%   1/h^k * [c(1)*u(x0 + j(1)*h) + ... + c(n)*u(x0 + j(n)*h)].
% Typically j(1) <= 0 <= j(n) and the values in j are monotonically
% increasing, but neither of these conditions is required.
% The routine fdcoeffF is used to compute the coefficients.
%
% Example:   fdstencil(2,-1:1);
%   determines the 2nd order centered approximation of the 2nd derivative.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

n = length(j);
if k>=n 
   error('***  length(j) must be larger than k')
end

c = fdcoeffF(k,0,j);      % coefficients for k'th derivative 

% print out stencil:

disp(' ')
fprintf('The derivative u^(%i) of u at x0 is approximated by',k)
disp(' ')
fprintf('1/h^%i * [',k)
for i=1:n-1
  if j(i) < 0
    fprintf('%22.15e * u(x0%i*h) + ',c(i),j(i))
  elseif j(i) == 0
    fprintf('%22.15e * u(x0) + ',c(i))
  else
    fprintf('%22.15e * u(x0+%i*h) + ',c(i),j(i))
  end
end

fprintf('%22.15e * u(x0+%i*h) ]   ',c(n),j(n))


% determine dominant terms in truncation error and print out:

err0 = c*(j(:).^n) / factorial(n);
err1 = c*(j(:).^(n+1)) / factorial(n+1);
if (abs(err0)) < 1e-14,  err0 = 0; end   % for centered approximations, expect
if (abs(err1)) < 1e-14,  err1 = 0; end   % one of these to be exactly 0.
disp(' ')
disp('For smooth u,')
fprintf('Error = %g * h^%i*u^(%i) + %g * h^%i*u^(%i) + ...',err0,n-k,n,err1,n-k+1,n+1)

disp(' ')