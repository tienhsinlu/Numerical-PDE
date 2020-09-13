% Euler Backward
T = 1;
f = @(u,t) -100*(u-cos(t));
[t,u_exact] = ode45(@(t,u) -100*(u-cos(t)), [0 1], 1);
dt = 0.1;
x = 0:dt:T;
u = zeros(1,length(x));
u(1)=1;
for i = 1:(length(x)-1)
    u(i+1)= (u(i)+ 100*cos(x(i)+dt)*dt)/(100*dt+1);
end
hold off;
plot(t,u_exact);
hold on;
scatter(x,u);