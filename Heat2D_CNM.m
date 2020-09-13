% 2D Heat Equation (Crank-Nicolson Method)

% parameters;
a_x=0; b_x=1;
a_y=0; b_y=1;
T=0.5;
u0=@(x,y)sin(pi*x).*sin(pi*y);
u_exact=@(x, y, t)exp(-2*pi^2*t).*sin(pi*x).*sin(pi*y);

M = [20,40,80,160,320];
% Storing data
time_CN = zeros(1,length(M));
error_CN_inf = zeros(1,length(M));
error_CN_2h = zeros(1,length(M));
for k = 1:length(M)
    m = M(k);
    % Grids;
    x=linspace(a_x, b_x, m+1);
    y=linspace(a_y, b_y, m+1);
    h=(b_x-a_x)/m;
    dt = h;
    [x_coord, y_coord]=meshgrid(x, y); 
    %x_coord, y_coord contains the x, y coordinates of each grid points.
    x_coord=x_coord';
    y_coord=y_coord';
    % Initialize
    U=u0(x_coord, y_coord);
    % Generating the semi-discrete problem
    I=eye(m-1);
    e=ones(m-1,1);
    mat_T=spdiags([e, -4*e, e], [-1,0,1], m-1, m-1);
    S=spdiags([e, e], [-1, 1], m-1, m-1);
    A=(kron(I, mat_T)+kron(S, I))/h^2;
    %Reordering the unknowns in the rowwise ordering and reshape the unknown matrix into a vector. 
    U = reshape(U(2:end-1,2:end-1),[],1);
    %Generate matrices I-dt/2*A and I+dt/2*A and store them in the sparse storage.
    LHS = speye((m-1)^2) - (dt/2)*A;
    RHS = speye((m-1)^2) + (dt/2)*A;
    
    %March in time and measure the execution time of the time marching.
    step = T/dt;
    tic;
    for i = 1:step
        Unew = LHS\(RHS*U);
        U = Unew;
    end
    time_CN(k) = toc();
    %Calculate the error.
    Exact = u_exact(x_coord,y_coord,T);
    Exact = reshape(Exact(2:end-1,2:end-1),[],1);
    error_CN_2h(k) = norm(U-Exact,2)*h;
    error_CN_inf(k) = norm(U-Exact,inf);
end
format short g
format compact
T = table(M',error_CN_2h',error_CN_inf',time_CN');
T.Properties.VariableNames = {'m','ErrorInf', 'Error2h','Time'};
disp("2D Heat Equation using Crank-Nicolson Method:");
disp(T);