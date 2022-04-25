clc; clear; close all;
%%  System Modeling & Simulation
%   2nd assignement, task 2, Lyapunov Mixed, 1 of 1 scripts

%%
a = 1.5;
b = 2;
u = @(t) 3*cos(2*t);
frequency = 30; 
amplitude = 0.25;
%frequencyRange = [0.5 1 2]; %[10 50 200]
%amplitudeRange =  [0.1 0.2 0.5 ]; %[1 2 5];
noise = @(t) amplitude*sin(2*pi*frequency*t);
gammaValues = [1 1];
thetaM = 0.5;
[tFromODE,x,xHat,aHat,bHat] = lyapunovMixed(u,gammaValues,thetaM,a,b, noise);
error = x - xHat;
plot(tFromODE, error)

% figure()
% subplot(2,3,1)
% plot(tFromODE,xHat,tFromODE,x),legend('$\hat{x}$','x','interpreter','latex')
% xlabel('Time (s)','interpreter','latex')
% title('$x \quad and \quad \hat{x}\quad values $','interpreter','latex')
% subplot(2,3,4)
% plot(tFromODE,x - xHat),legend('$\hat{x} - x$','interpreter','latex')
% xlabel('Time (s)','interpreter','latex')
% title('$x - \hat{x} \quad error $','interpreter','latex');
% 
% subplot(2,3,2)
% plot(tFromODE,aHat,tFromODE, a * ones(1, length (tFromODE))),legend('$\hat{a}$','a','interpreter','latex')
% xlabel('Time (s)','interpreter','latex')
% title('$a \quad and \quad \hat{a}\quad values $','interpreter','latex')
% subplot(2,3,5)
% plot(tFromODE, a - aHat),legend('$\hat{a} - a$','interpreter','latex'),xlabel('Time (s)','interpreter','latex')
% title('$\hat{\alpha}$','interpreter','latex')
% title('$a - \hat{a} \quad error $','interpreter','latex');
% 
% subplot(2,3,3)
% plot(tFromODE,bHat,tFromODE, b * ones(1, length (tFromODE))),legend('$\hat{a}$','a','interpreter','latex')
% xlabel('Time (s)')
% title('$b \quad and \quad \hat{b}\quad values $','interpreter','latex')
% subplot(2,3,6)
% plot(tFromODE, b - bHat),legend('$b - \hat{b} $','interpreter','latex'),xlabel('Time (s)','interpreter','latex')
% title('$b - \hat{b} \quad error $','interpreter','latex');
% 

function [tFromODE,x,xHat,aHat,bHat] = lyapunovMixed(u,gamma1gamma2,thetaM,a,b,noise)
equationsSystem = @(t,y) prepareEquationsSystemToBeSolved_MixedLyapunov(t,u,y,gamma1gamma2,thetaM,a,b,noise);
[tFromODE,y] = ode45(equationsSystem, [0 100], [0; 0; 0; 0]);
x = y(:,1);
xHat = y(:,2);
aHat = y(:,3);
bHat = y(:,4);
end

function equations = prepareEquationsSystemToBeSolved_MixedLyapunov(t,u,y,gamma1gamma2,thetaM,a,b,noise)
x = y(1); xHat=y(2); theta1 = y(3); theta2 = y(4);
xDot = -a*x + b*u(t);
error = x + noise(t) - xHat;
xHatDot = -theta1*(x + noise(t)) + theta2*u(t) + thetaM*error;
thetaDot1 = -gamma1gamma2(1)*error*(x + noise(t));
thetaDot2 = gamma1gamma2(2)*error*u(t);  

equations = [xDot; xHatDot; thetaDot1; thetaDot2];

end

