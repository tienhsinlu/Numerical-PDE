%spatial domain [0, 2*pi]
a=0;
b=2*pi;

%initial condition (periodic function)
u0=@(x) sin(x);

%time step size
dt=0.05;
n_t=500;
T=n_t*dt;

figure;
%A uniform grid with 101 points are placed over the spatial domain.
for j=1:6
    m=2^j*10;
    %Generate the grid
    x=linspace(a,b,m+1);
    
    %grid/mesh size
    dx=(b-a)/m;
    
    
    u=u0(x);
    for i=2:n_t
       u(2:m)=0.5*(u(3:m+1)+u(1:m-1))-dt/(2*dx)*(u(3:m+1)-u(1:m-1));
       u(1)=0.5*(u(2)+u(end))-dt/(2*dx)*(u(2)-u(end));
       u(end)=0.5*(u(end-1)+u(1))-dt/(2*dx)*(u(1)-u(end-1));
    end
    
    u_exact=u0(x-T);
    
    subplot(2,3,j)
    plot(x, u, 'bo', 'LineWidth', 2);
    hold on
    plot(x, u_exact, 'r--', 'LineWidth', 2);
    hold off
    xlabel('x');
    ylabel('u');
    legend('Approximation', 'Reference', 'Location', 'best');
    title(['m=', num2str(m), ', lambda=', num2str(dt/dx)])
end