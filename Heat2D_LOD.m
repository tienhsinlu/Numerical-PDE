% 2D Heat Equation (LOD)

% parameters;
a_x=0; b_x=1;
a_y=0; b_y=1;
T=0.5;
u0=@(x,y)sin(pi*x).*sin(pi*y);
u_exact=@(x, y, t)exp(-2*pi^2*t).*sin(pi*x).*sin(pi*y);

M = [20,40,80,160,320];
% Storing data
time_LOD = zeros(1,length(M));
error_LOD_2h = zeros(1,length(M));
error_LOD_inf = zeros(1,length(M));
for k = 1:length(M)
    m = M(k);
    % Grids;
    x=linspace(a_x, b_x, m+1);
    y=linspace(a_y, b_y, m+1);
    h=(b_x-a_x)/m;
    dt = T/m;
    [x_coord, y_coord]=meshgrid(x, y); 
    %x_coord, y_coord contains the x, y coordinates of each grid points.
    x_coord=x_coord';
    y_coord=y_coord';
    % Initialize
    U=u0(x_coord, y_coord);
    % Generating the semi-discrete problem
    A = diag(-2*ones(1,m-1))+diag(ones(1,m-2),1)+diag(ones(1,m-2),-1);
    A = (1/h^2)*A;
    A = sparse(A);
    LHS = speye(m-1)-(dt/2)*A;
    RHS = speye(m-1)+(dt/2)*A;
    U = U(2:end-1,2:end-1);
    U = sparse(U);
    %March in time and measure the execution time.
    step = T/dt;
    U_temp = zeros(m-1,m-1);
    U_new = zeros(m-1,m-1);
    tic;
    for i=1:step
        U_temp = LHS\(RHS*U);
        U_new = LHS\(RHS*U_temp);
        U = U_new;
    end
    time_LOD(k) = toc();
    %Calculate the error.
    Exact = u_exact(x_coord,y_coord,T);
    U = reshape(U,[],1);
    Exact = reshape(Exact(2:end-1,2:end-1),[],1);
    error_LOD_2h(k) = norm(U-Exact,2)*h;
    error_LOD_inf(k) = norm(U-Exact,inf);
end
format short g
format compact
T = table(M',error_LOD_2h',error_LOD_inf',time_LOD');
T.Properties.VariableNames = {'m','ErrorInf', 'Error2h','Time'};
disp("2D Heat Equation using LOD Method:");
disp(T);