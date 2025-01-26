clear
%% 4 NMDE
nx = 100;
G = numgrid('S',nx);
A = (delsq(G)*(nx-1)^2);
n = length(A);
lambda = -eigs(A,1,'lm');

hbar_rk = -2.78/(2*lambda); 
fprintf('hbar_rk %9.4e' ,hbar_rk);

%% obtain y0 and y_exact:
y0 = ones(n,1);
yex = load('exact_solution.txt','r');

%% Computing ode45 Method:
tspan = [0 0.1];

tStart = tic;  
[t,y] = ode45(@(t,y) func_A(t,y,A), tspan, y0);
tEnd = toc(tStart);

CPU_time_ode45 = tEnd; 
num_steps_ode45 = length(t);

%% Computation of the error:
err_ode45 = max(abs(y(end,:)'- yex(:,2)));

%% display(err_ode45);
 run1=["run 1", "ode45", num_steps_ode45, err_ode45, CPU_time_ode45];
 display(run1);

%% Computing Crank Nicolson Method:
hh = [0.001 0.0001 0.00001];
err_cn = zeros(3,1);
num_steps_cn = zeros(3,1);
CPU_time_cn = zeros(3,1);

for j = 1:3
   my_h = hh(j);
   my_A = (speye(n)+0.5*my_h*A);
   num_times = round(0.1/my_h);
   y = zeros(n,num_times+1);
   y(:,1)=y0;
   tStart = tic;
    for i=2:num_times+1
        b = (speye(n)-0.5*my_h*A)*y(:,i-1);
        y(:,i)= pcg(my_A, b, my_h^3,100000);
    end
   tEnd = toc(tStart);
   num_steps_cn(j) = num_times;
   CPU_time_cn(j) = tEnd;

%% Computation of the error
err_cn(j) = max(abs(y(:,end)-yex(:,2)));

ycn_due(:,j) = y(:,2);
ycn_tre(:,j) = y(:,3);

end

run2=["run 2", "CN method", num_steps_cn(1), err_cn(1), CPU_time_cn(1)];
run3=["run 3", "CN method", num_steps_cn(2), err_cn(2), CPU_time_cn(2)];
run4=["run 4", "CN method", num_steps_cn(3), err_cn(3), CPU_time_cn(3)];

display(run2);
display(run3);
display(run4);

%% Computing BDF3 Method:
hh = [0.001 0.0001 0.00001];
err_BDF3 = zeros(3, 1);
num_steps_BDF3 = zeros(3, 1);
CPU_time_BDF3 = zeros(3, 1);

for j=1:3
my_h = hh(j);
num_times = round(0.1/my_h);
my_A = (speye(n)+6/11*my_h*A);

y = zeros(n,num_times+1);
y(:,1)=y0;
y(:,2)=ycn_due(:,j);
y(:,3)=ycn_tre(:,j);
tStart = tic;
for i=4:num_times+1
    b = +18/11*y(:,i-1)-9/11*y(:,i-2)+2/11*y(:,i-3);
    y(:,i)= pcg(my_A,b, my_h^3, 100000);
end
tEnd = toc(tStart);
CPU_time_BDF3(j) = tEnd;
num_steps_BDF3(j) = num_times;

%% Computation of the error
err_BDF3(j) = max(abs(y(:,end)-yex(:,2)));
end

run5=["run 5", "BDF3", num_steps_BDF3(1), err_BDF3(1), CPU_time_BDF3(1)];
run6=["run 6", "BDF3", num_steps_BDF3(2), err_BDF3(2), CPU_time_BDF3(2)];
run7=["run 7", "BDF3", num_steps_BDF3(3), err_BDF3(3), CPU_time_BDF3(3)];

display(run5);
display(run6);
display(run7);

%% Final table
head=["RUN", "METHOD", "NUMBER OF STEPS", "ERROR", "CPU TIME"];
TABLE = [head; run1; run2; run3; run4; run5; run6; run7];
display(TABLE);


%% Function for ode45
function dydt = func_A(t,y,A)
dydt = -A*y;
end
