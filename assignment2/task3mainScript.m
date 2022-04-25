clc; clear; close all;
%%  System Modeling & Simulation
%   2nd assignement, task 3, main script, 1 of 1 scripts

%% MIMO
% the parameters I want to estimate
a = [-0.5 -3; 4 -2];
b = [1;1.4];
% input
u = @(t) 7.5*cos(3*t) + 10*cos(2*t);
gammaRange = 1; 
gamma = [1 1];
[tFromODE,x,xHat,aHat,bHat] = lyapunovParallelMIMO(u,a,b,gamma);

%% plots below
a = a';
figure
for j=1:4
    hold on
    yline(a(j),'r--','HandleVisibility','off');
    plot(tFromODE,aHat(:,j,1))
    legend('$\hat{a}_{1,1}$','$\hat{a}_{1,2}$','$\hat{a}_{2,1}$','$\hat{a}_{2,2}$','interpreter','latex')
    xlabel('Time (s)')
end
figure
for j=1:2
    hold on
    yline(b(j),'r--','HandleVisibility','off');
    plot(tFromODE,bHat(:,j,1))
    xlabel('Time (s)')
    legend('$\hat{b}_{1}$','$\hat{b}_{2}$','interpreter','latex')
end

figure
subplot(1,2,1)
plot(tFromODE,xHat(:,1,1),tFromODE,x(:,1))
legend('$\hat{x_1}$', '$x_1$', 'interpreter', 'latex')
xlabel('Time (s)')
subplot(1,2,2),plot(tFromODE,xHat(:,2,1),tFromODE,x(:,2))
legend('$\hat{x_2}$', '$x_2$', 'interpreter', 'latex')
xlabel('Time (s)')

function [tFromODE,x,xHat,AHat,BHat] = lyapunovParallelMIMO(u,a,b,gamma)
equationsSystem = @(t,y) prepareEquationsSystemToBeSolved_Parallel_MIMO_Lyapunov(t,u,y,gamma,a,b);
[tFromODE,y] = ode45(equationsSystem, [0 100], zeros(10,1));
x    = [y(:,1) y(:,2)];
xHat = [y(:,3) y(:,4)];
AHat = [y(:,5) y(:,6) y(:,7) y(:,8)];
BHat = [y(:,9) y(:,10)];
end

function equations = prepareEquationsSystemToBeSolved_Parallel_MIMO_Lyapunov(t,u,y,gamma,a,b)
x    = [y(1); y(2)];
xHat = [y(3); y(4)];
aHat = [y(5) y(6); y(7) y(8)];
bHat = [y(9); y(10)];

xDot = a*x + b*u(t);

error = [x(1) - xHat(1); x(2) - xHat(2)];
xHatDot = aHat*xHat + bHat*u(t);
aHatDot = gamma(1) * (error*xHat');
bHatDot = gamma(2)*error*u(t);

equations = [xDot(1); xDot(2); xHatDot(1); xHatDot(2);...
    aHatDot(1); aHatDot(3); aHatDot(2); aHatDot(4);...
    bHatDot(1); bHatDot(2)];
end