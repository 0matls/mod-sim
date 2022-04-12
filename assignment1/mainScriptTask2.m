clear; clc; close all;
format long;

%%  System Modeling & Simulation
%   1st assignement, task 2, main script + 1 function at the bottom

%%
% sampling based on the provided v.p 
samplingTimes = 0:0.000001:5;
samples = length(samplingTimes);
VrSamples = NaN(1,samples);
VcSamples = NaN(1,samples);

for i = 1:samples
    [VrSamples(1,i), VcSamples(1,i)] = v(samplingTimes(i));
end

% the 2 inputs
input1 = @(t)3*sin(2*t);
input2 = ones(1,samples);
input2 = input2 + 1;

%% initial plots
% figure();
% plot(samplingTimes,VcSamples, 'color','magenta');
% title("V_C");
% xlabel("Time");
% ylabel("Volt");
% 
% figure();
% plot(samplingTimes,VrSamples, 'color', 'magenta');
% title("V_R");
% xlabel("Time");
% ylabel("Volt");

%%
polesRange = [200]; %[25 40 50 75 100 200 500 1000]; % testing poles
poleIndex = length(polesRange);
theta = NaN(6,poleIndex);
lambdaVector = NaN(2,poleIndex);
count = 0;

for i =1 : poleIndex
    count = count + 1;
    fprintf("%d of %d\n",count, poleIndex);
    pole1 = -polesRange(count);
    pole2 = -polesRange(count);
    lambdaVector(:,count) = [-(pole1+pole2) ; pole1*pole2];
    theta(:,count) = thetaEstimator2(pole1, pole2, samples, VcSamples, input1, input2, samplingTimes);
    RC1Inversed(count) = theta(1,count) + lambdaVector(1,count);
    LC1Inversed(count) = theta(2,count) + lambdaVector(2,count);
    RC2Inversed(count) = theta(3,count);
    RC3Inversed(count) = theta(5,count);
    LC2Inversed(count) = theta(6,count);
    
    % calculate initial conditions
    dV = v(1e-10);
    Vzero = v(0);
    x0(1) = Vzero(1);
    x0(2) = (dV(1) - Vzero(1))/1e-10; % approximate derivative

    % ode values
    [t,x] = ode45(@(time,x)rlcODE(time,x, RC1Inversed(end), LC1Inversed(end)), samplingTimes, x0);

    estimates1(count) = abs(VcSamples(end) - x(end,1));
    VRhat = 2 + input1(t(end))' - x(end,1)';
    estimates2(count) = abs(VRhat);

end


%% plots to support poles' decision
% figure();
% plot(polesRange, estimates1);
% title("Estimation Error of V_C");
% grid on;
% figure();
% plot(polesRange, estimates2);
% grid on;
% title("Estimation Error of V_R");

% figure();
% plot(polesRange, RC1Inversed, 'color', 'green');
% title("1/RC parameters convergence ");
% hold on
% plot(polesRange, RC2Inversed,  'color', 'blue');
% plot(polesRange, RC3Inversed, 'color', 'magenta');
% xlabel('poles'' values of $\Lambda(s)$', 'interpreter', 'latex');
% 
% figure();
% hold on
% plot(polesRange, LC1Inversed, 'color', 'magenta'); 
% title("1/LC parameters convergence");
% plot(polesRange, LC2Inversed, 'color', 'green');
% xlabel('poles''s values of $\Lambda(s)$', 'interpreter', 'latex');


%% error VC - when I have just one (double) pole, not multiple in a vector
% figure();
% plot(samplingTimes, VcSamples(:) - x(:,1));
% title("Estimation Error of VC at p = -200");


%% part 2 below - add three outliers 
outlierIndex = [round(samples/60) round(samples/2) round(samples/1.1)];
outlier = [10000 4500 35000]; % or 100 100 100
VCnoised = VcSamples;
VCnoised(outlierIndex) = VCnoised(outlierIndex) + outlier;
VRnoised = VrSamples;
VRnoised(outlierIndex) = VRnoised(outlierIndex) + outlier;

figure();
plot(samplingTimes,VCnoised);
title("V_C with noised added because of outliers");
xlabel("Time");
ylabel("Volt");

% figure();
% plot(samplingTimes,VCnoised);
% title("V_R with noised added because of outliers");
% xlabel("Time");
% ylabel("Volt");

theta = thetaEstimator2(-50, -50, samples, VCnoised, input1, input2, samplingTimes);

invRC1 = theta(1) + 100
invLC1 = theta(2) + 2500
invRC2 = theta(3);
invRC3 = theta(5);
invLC2 = theta(6);

fprintf('1/RC1=%f\n1/RC2=%f\n1/RC3=%f\n1/LC1=%f\n1/LC2=%f\n', invRC1,... 
        invRC2, invRC3, invLC1, invLC2);
   
% calculate initial conditions
dV = v(1e-10);
Vzero = v(0);
x0(1) = Vzero(1);
x0(2) = (dV(1) - Vzero(1))/1e-10; % approximated derivative

[t,x] = ode45(@(time,x)rlcODE(time, x, invRC1, invLC1), samplingTimes, x0);

% error 
figure();
plot(samplingTimes, VcSamples(:) - x(:,1));
title("Estimation error of V_C with 3 outliers");


function ds = rlcODE(t,x,invRC,invLC)

ds(1) = x(2);
ds(2) = -x(2)*invRC - x(1)*invLC + 6*cos(2*t)*invRC + 0 + 2*invLC;
ds = ds';

end

