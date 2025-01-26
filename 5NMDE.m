clear
%% 5 NMDE
t0 = 0; 
T = 300;

%% Timesteps:
h = 0.001;
N = (T-t0)/h;

%% Value of parameters:
alpha=0.2;
beta=0.01;
gamma=0.004;
delta=0.07;

%% Functions:
f = @(x, y)  x*(alpha-(beta*y));
g = @(x, y)  y*((gamma*x)-delta);

%%  computing Runge Kutta method:
xin(1)=19;
yin(1)=22;

for n = 1:N

    K1x = f(xin(n), yin(n));
    K1y = g(xin(n), yin(n));
    K2x = f(xin(n) + h*K1x/2, yin(n) + h*K1y/2);
    K2y = g(xin(n) + h*K1x/2, yin(n) + h*K1y/2);
    K3x = f(xin(n) + h*K2x/2, yin(n) + h*K2y/2);
    K3y = g(xin(n) + h*K2x/2, yin(n) + h*K2y/2);
    K4x = f(xin(n) + h*K3x, yin(n) + h*K3y);
    K4y = g(xin(n) + h*K3x, yin(n) + h*K3y);

    xin(n+1) = xin(n) + h*(K1x + 2*K2x + 2*K3x + K4x)/6;
    yin(n+1) = yin(n) + h*(K1y + 2*K2y + 2*K3y + K4y)/6;

end

%% Plotting the value:
tt= t0:h:T;

 figure(1)
 hold on
 plot(tt,xin,'*r-')
 plot(tt,yin,'*c-')
 legend('Prey', 'Predator')
 hold off 
