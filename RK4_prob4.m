% Using RK-4 to solve the ODE with dt=0.1, and 0.001;
T = 1;
f = @(u,t) -100*(u-cos(t));
[t,u_exact] = ode45(@(t,u) -100*(u-cos(t)), [0 1], 1);
dt = 0.1;
x = 0:dt:T;
u = zeros(1,length(x));
u(1)=1;
for j = 1:(length(x)-1)
    tn = x(j);
    Y1 = u(j);
    Y2 = Y1 + dt/2*f(Y1,tn);
    Y3 = Y1 + dt/2*f(Y2,tn+dt/2);
    Y4 = Y1 + dt*f(Y3,tn+dt/2);

    u(j+1) = u(j) + dt/6*(f(Y1,tn)+2*f(Y2,tn+dt/2)+2*f(Y3,tn+dt/2)+...
        f(Y4,tn+dt));
end
hold off;
plot(t,u_exact);
hold on;
scatter(x,u);