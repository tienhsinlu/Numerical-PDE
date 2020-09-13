% Homework 9 Problem 4: Implement the scheme
a_x=0; b_x=1;
a_y=0; b_y=1;
f = @(x,y) (2*pi^2+1)*sin(pi.*x).*sin(pi.*y);
u_exact = @(x,y) sin(pi.*x).*sin(pi.*y);
M = [80,160,320,640];
error = zeros(4,1);
for k = 1:length(M)
    m = M(k);
    x=linspace(a_x, b_x, m+1);
    y=linspace(a_y, b_y, m+1);
    h=(b_x-a_x)/m;
    [x_coord, y_coord]=meshgrid(x, y); 
    x_coord=x_coord';
    y_coord=y_coord';
    % Generate the matrix
    I=eye(m-1);
    e=ones(m-1,1);
    mat_T=spdiags([1*e, -(4+h^2)*e, 1*e], [-1,0,1], m-1, m-1);
    S=spdiags([e, e], [-1, 1], m-1, m-1);
    A = (kron(I, mat_T)+kron(S, I));
    A = -1/h^2*A;
    b = f(x_coord,y_coord);
    b = reshape(b(2:end-1,2:end-1),[],1);
    % Use CG
    U = cg(A,b,10e-6,1000); % function cg is attached
    Exact = u_exact(x_coord,y_coord);
    Exact = reshape(Exact(2:end-1,2:end-1),[],1);
    error(k) = norm(U-Exact,2)*h;
end
format short g
format compact
T = table(M',error);
T.Properties.VariableNames = {'m','Error'};
disp("Error:");
disp(T);