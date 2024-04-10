clc; clear; close all;
%%  System Modeling & Simulation
% project - main script - 2 of 2 scripts 
%%
a = 2.3;
b = 1.2;
xDesired = @(t) 0.5*cos(3*t);  
w = 2;
rZero = 4;
rInf = 0.01; %fixed 
lambda = 2; %fixed
k = 5;
[tFromODE,x,xMine] = lyapunovMixed(xDesired,a,b,w,rZero,rInf,lambda,k);
%%
rho = @(t)  (rZero - rInf)*exp(-lambda*t) + rInf;
z = @(t,x) (x - xDesired(t)) ./ rho(t);
epsilonFunction = @(t,x) log( (1 + z(t,x)) / (1 - z(t,x)) ); 
r = @(t,x) 2 ./ rho(t) ./ (1 - z(t,x).^2 );
u = @(t,w,x) -k./w .* (r(t,x) + 1./r(t,x)).*epsilonFunction(t,x);
error = x - xMine;
error_x_x_d = x - xDesired(tFromODE);
%% plots
figure(1)
subplot(3,1,1)
plot(tFromODE, error_x_x_d)
xlabel('Time (s)','interpreter','latex')
title('$-\rho(t) < x(t) -x_{d}(t) < \rho(t) $','interpreter','latex');
hold on
plot(tFromODE, rho(tFromODE))
plot(tFromODE, -rho(tFromODE))
legend("x - x_{desired}", "rho(t)", "-rho(t)") 
subplot(3,1,2)
plot(tFromODE, x);
xlabel('Time (s)','interpreter','latex')
hold on 
plot(tFromODE, xMine);
plot(tFromODE, error);
legend("x", "xHat", "error") 
title('$x - \hat{x} \quad and \quad their \quad error $','interpreter','latex');
subplot(3,1,3)
plot(tFromODE, u(tFromODE,w,x))
xlabel('Time (s)','interpreter','latex')
title('$u(t) $','interpreter','latex');
%%
function [tFromODE,x,xMine] = lyapunovMixed(xDesired,a,b,w,rZero,rInf,lambda,k)
rho = @(t)  (rZero - rInf)*exp(-lambda*t) + rInf;
z = @(t,x) (x - xDesired(t))/rho(t);
epsilonFunction = @(t,x) log( (1 + z(t,x)) / (1 - z(t,x)) ); 
r = @(t,x) 2 / rho(t) / (1 - z(t,x)^2 );
u = @(t,x) -k/w * (r(t,x) + 1/r(t,x))*epsilonFunction(t,x);
equationsSystem = @(t,y) prepareEquationsSystemToBeSolved_MixedLyapunov(t,y,u,a,b);
[tFromODE,y] = ode15s(equationsSystem, [0 30], [0; 0]);
x = y(:,1);
xMine = y(:,2);
end

function equations = prepareEquationsSystemToBeSolved_MixedLyapunov(t,y,u,a,b)
x = y(1); xMine=y(2); 
xDot = sys(x, u(t,x));
xDotMine = a*xMine + b*u(t,xMine);
equations = [xDot; xDotMine];
end
