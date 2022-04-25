clc; clear; close all;

%%  System Modeling & Simulation
%   2nd assignement, task 2, Lyapunov Parallel, 1 of 1 scripts

a = 1.5;
b = 2;
u = @(t) 3*cos(2*t);
frequency = 30; 
amplitude = 0.25;
noise = @(t) amplitude*sin(2*pi*frequency*t);
%gammaRange = [0.5 1 5 10 30];
gammaValues = [5 5];
%frequencyRange =  [10 50 200]; %[0.5 1 2]
% amplitudeRange =  [1 2 5]; %[0.1 0.2 0.5 ];

[tFromODE,x,xHat,aHat,bHat] = lyapunovParallel(u,gammaValues,a,b, noise);
error = x - xHat;

figure
plot(tFromODE, error)

function [tFromODE,x,xHat,aHat,bHat] = lyapunovParallel(u,gamma1gamma2,a,b,noise)
equationsSystem = @(t,y) prepareEquationsSystemToBeSolved_ParallelLyapunov(t,u,y,gamma1gamma2,a,b, noise);
[tFromODE,y] = ode45(equationsSystem, [0 80], [0; 0; 0; 0]);
x = y(:,1);
xHat = y(:,2);
aHat = y(:,3);
bHat = y(:,4);
end

function equations = prepareEquationsSystemToBeSolved_ParallelLyapunov(t,u,y,gamma1gamma2,a,b, noise)
x = y(1); xHat=y(2); theta1 = y(3); theta2 = y(4);
xDot = -a*x + b*u(t);
xHatDot = -theta1*xHat + theta2*u(t);
error = x + noise(t) - xHat;
thetaDot1 = -gamma1gamma2(1)*error*xHat;
thetaDot2 = gamma1gamma2(2)*error*u(t);

equations = [xDot; xHatDot; thetaDot1; thetaDot2];
end
