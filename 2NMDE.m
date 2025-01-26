clear
%% 2 NMDE
lambda = 10;
t0 = 0;
T = 2;
%% equation:
%% y' = -10*y^2
%% y(0) = 1
f = @(y)-lambda*(y.^2);
y(1)= 1;
errs_fin= ones(6,1);

%% exact solution:
true_sol = @(t) 1./(1+lambda*t);

for M=5:1:10

h= 2^(-M);
N= (T-t0)/h;
tt= t0:h:T;
yy= true_sol(tt);
          
%% compute the 4th order Runge-Kutta method:
for i= 1:1:N

     yrk(1)= y(1);
     K1= f(yrk(i));
     K2= f(yrk(i) + h*K1/2);
     K3= f(yrk(i) + h*K2/2);
     K4= f(yrk(i) + h*K3);
     yrk(i+1)= yrk(i) + h*(K1 + 2*K2 + 2*K3 + K4)/6;

end

%% Computing the final errors:

 errs_fin(M-4)= abs(yy(N+1)-yrk(N+1));

end

%% plotting the final errors:

figure(1)
axes('Xscale', 'log', 'Yscale', 'log')
xlabel('Number of steps N');
ylabel('Final errors');
hold on

for M= 5:1:10

    h= 2^(-M);
    N= (T-t0)/h;
    plot(N',errs_fin(M-4)','*r-','LineWidth',1)
end

hold off

