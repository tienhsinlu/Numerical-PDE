% Order of Accuracy Analysis
T=0.5;
lambda=0.6;
n_refine=6;
err_l2dx=zeros(n_refine,1);

order=@(e) log(e(1:end-1)./e(2:end))./log(2);
%spatial domain [0, 2*pi]
a=0;
b=2*pi;
%initial condition (periodic function)
u0=@(x) sin(x);

for j=1:n_refine
    m=2^j*10;
    x=linspace(a,b,m+1);
    dx=(b-a)/m;
    dt=lambda*dx;
    n_t=floor(T/dt);
    
    u=u0(x);
    for i=2:n_t
       u(2:m)=0.5*(u(3:m+1)+u(1:m-1))-dt/(2*dx)*(u(3:m+1)-u(1:m-1));
       u(1)=0.5*(u(2)+u(end))-dt/(2*dx)*(u(2)-u(end));
       u(end)=0.5*(u(end-1)+u(1))-dt/(2*dx)*(u(1)-u(end-1));

    end
    
    %after the n_t steps, we may still have T-n_t*dt to arrive at the final time T.
    dt=T-n_t*dt;
    u(2:m)=0.5*(u(3:m+1)+u(1:m-1))-dt/(2*dx)*(u(3:m+1)-u(1:m-1));
    u(1)=0.5*(u(2)+u(end))-dt/(2*dx)*(u(2)-u(end));
    u(end)=0.5*(u(end-1)+u(1))-dt/(2*dx)*(u(1)-u(end-1));
    u_exact=u0(x-T);
    
    err_l2dx(j)=norm(u-u_exact, 2)*sqrt(dx);
end

format short e
disp('Numerical errors and orders in the l_{2,dx} norm:');
disp([err_l2dx, [0; order(err_l2dx)]])
