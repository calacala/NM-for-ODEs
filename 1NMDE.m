clear
%% 1 NMDE
lambda = 5;
t0 = 0; 
T = 10;

%% equation: 
%% y' = -5*y 
%% y(0) = 1
f = @(y)-lambda*y;

%% Timesteps:
h = 0.02;
N = (T-t0)/h; 
T = h*N+t0;

%% exact solution:
true_sol = @(t) exp(-5*t);
tt = t0:h:T;
yy = true_sol(tt);

%% Simpson's Method, obtained y(n+1) with Farward Euler Method:
yfe(1) = 1;
yfe(2) = yfe(1) + h*f(yfe(1));

for n = 1:N
    yfe(n+2) = fzero(@(x) x - yfe(n) - (h/3)*(f(yfe(n))+4*f(yfe(n+1))+f(x)), yfe(n)); 
end

%% Simpson's Method, obtained y(n+1) with 4th order Runge Kutta Method:
yrk(1) = 1;

    K1 = f(yrk(1));
    K2 = f(yrk(1) + h*K1/2);
    K3 = f(yrk(1) + h*K2/2);
    K4 = f(yrk(1) + h*K3);
    yrk(2) = yrk(1) + h*(K1 + 2*K2 + 2*K3 + K4)/6;

for n = 1:N
    yrk(n+2) = fzero(@(z)  z - yrk(n) - (h/3)*(f(yrk(n))+4*f(yrk(n+1))+f(z)), yrk(n)); 
end

%% plotting the errors:
yfe(end) = [];
yrk(end) = [];
err_fe = abs(yfe-yy);
err_rk = abs(yrk-yy);
display(err_fe);
display(err_rk);

figure(1)
hold on
plot(tt,err_fe,'*r-', 'Linewidth', 1)
tvec = linspace(t0,T);
legend('Forward Euler Method')
hold off 

figure(2)
hold on
plot(tt,err_rk,'*c-', 'Linewidth', 1)
tvec = linspace(t0,T);
legend('4th Runge Kutta Method')
hold off 