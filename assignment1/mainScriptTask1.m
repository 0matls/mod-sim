clear; clc; close all;

%%  System Modeling & Simulation
%   1st assignement, task 1, main script

%% given parameters and values
m = 10;
b = 0.3;
k = 1.5;
u = @(t) 10*sin(3*t) + 5;

%% time-sampling and initial contitions
tspan = 0:0.1:10;
state0 = [0 0]; % zero initial conditions

%%
f = @(t,s)msdODE(s, t, k, b, m, u);
options = odeset('AbsTol', 10^(-11), 'RelTol', 10^(-10));
[time, y] = ode45(f, tspan, state0, options);

%% poles optimization
rangeStart = 1;
rangeEnd = 8;
rangeStep = 0.1;
polesRange = rangeStart:rangeStep:rangeEnd; 

jIndex = 0;
iIndex = 0;
polesIndex = 0;

for i = polesRange
    iIndex = iIndex  + 1;
    for j = polesRange
    polesIndex = polesIndex + 1;
    jIndex = jIndex + 1;
    pole1 = -i;
    pole2 = -j;
    lambdaVector(:,polesIndex) = [-(pole1+pole2) pole1*pole2]; % the filter's coefficients
    thetaEstimation(:,polesIndex) = thetaEstimator(pole1, pole2, y, u, tspan);
    
    %% find the theta vector
    theta(1,iIndex,jIndex) = thetaEstimation(1,polesIndex) + lambdaVector(1,polesIndex);
    theta(2,iIndex,jIndex) = thetaEstimation(2,polesIndex) + lambdaVector(2,polesIndex);
    theta(3,iIndex,jIndex) = thetaEstimation(3,polesIndex);
    %% find the estimetions of the parameters
    mEstimated(iIndex,jIndex) = 1/theta(3,iIndex,jIndex);
    kEstimated(iIndex,jIndex) = theta(2,iIndex,jIndex) * mEstimated(iIndex,jIndex);
    bEstimated(iIndex,jIndex) = theta(1,iIndex,jIndex) * mEstimated(iIndex,jIndex);
    end
    jIndex = 0;
end


figure();
title("Estimation error altogether");
Z =  -abs(m - mEstimated) - abs(k - kEstimated) - abs(b - bEstimated);
surf(polesRange,polesRange,Z);
grid on;
legend("sum of errors");
xlabel('poles $\Lambda(s)$', 'interpreter','latex');
[Zmax,Idx] = max(Z(:))
[ZmaxRow,ZmaxCol] = ind2sub(size(Z), Idx)
pole1optimized = rangeStart + ZmaxRow*rangeStep
pole2optimized = rangeStart + ZmaxCol*rangeStep
mE = mEstimated(ZmaxRow,ZmaxCol)
kE = kEstimated(ZmaxRow,ZmaxCol)
bE = bEstimated(ZmaxRow,ZmaxCol)
