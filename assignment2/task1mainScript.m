clc; clear; close all;

%%  System Modeling & Simulation
%   2nd assignement, task 1, main script, 1 of 1 scripts

%% gradient descend method implementation for parameter estimmation
% the 2 different inputs u = 3 and u = 3cos(2t)
% u = @(t) 3*cos(0*t); % basically u = 3 
u = @(t) 3*cos(0*t);

% parameters a and b to be estimated
a = 1.5;
b = 2;

% the pole of the Ë(s) filter
poleRange = -1.9;
%poleRange = [-0.5 -1 -1.5 -2 -3];
%gammaRange = [0.5 1 5 30 200 10000];
gamma = 34;
% counter = 1;
% figure()
% hold on
for i=poleRange
    [tspan, x, phi1, phi2,aHat,bHat] = gradientDescent(u,i,gamma,a,b);
    theta1 =  -i -aHat;
    theta2 = bHat;
    xHat = theta1.*phi1 + theta2.*phi2;
     
%     legendInfo{counter} = ['pole = ' num2str(i)]; 
%     plot(tspan,x-xHat),legend(legendInfo,'interpreter','latex'),xlabel('Time (s)','interpreter','latex')
%     yline(0,'r--','HandleVisibility','off');
%     title('$x - \hat{x} \quad error $','interpreter','latex');
%     plot(tspan,a-aHat),legend(legendInfo1,'interpreter','latex'),xlabel('Time (s)','interpreter','latex')
%     yline(0,'r--','HandleVisibility','off');
%     title('$a - \hat{\alpha} \quad error $','interpreter','latex');
%     plot(tspan,b-bHat),legend(legendInfo2,'interpreter','latex'),xlabel('Time (s)','interpreter','latex')
%     yline(0,'r--','HandleVisibility','off');
%     title('$b - \hat{b}\quad error $','interpreter','latex');
%     counter = counter + 1;
end

figure()
subplot(2,3,1)
plot(tspan,xHat,tspan,x),legend('$\hat{x}$','x','interpreter','latex')
xlabel('Time (s)','interpreter','latex')
title('$x \quad and \quad \hat{x}\quad values $','interpreter','latex')
subplot(2,3,4)
plot(tspan,x - xHat),legend('$\hat{x} - x$','interpreter','latex')
xlabel('Time (s)','interpreter','latex')
title('$x - \hat{x} \quad error $','interpreter','latex');

subplot(2,3,2)
plot(tspan,aHat,tspan, a * ones(1, length (tspan))),legend('$\hat{a}$','a','interpreter','latex')
xlabel('Time (s)','interpreter','latex')
title('$a \quad and \quad \hat{a}\quad values $','interpreter','latex')
subplot(2,3,5)
plot(tspan, a - aHat),legend('$\hat{a} - a$','interpreter','latex'),xlabel('Time (s)','interpreter','latex')
title('$\hat{\alpha}$','interpreter','latex')
title('$a - \hat{a} \quad error $','interpreter','latex');

subplot(2,3,3)
plot(tspan,bHat,tspan, b * ones(1, length (tspan))),legend('$\hat{a}$','a','interpreter','latex')
xlabel('Time (s)')
title('$b \quad and \quad \hat{b}\quad values $','interpreter','latex')
subplot(2,3,6)
plot(tspan, b - bHat),legend('$b - \hat{b} $','interpreter','latex'),xlabel('Time (s)','interpreter','latex')
title('$b - \hat{b} \quad error $','interpreter','latex');

function [tFromODE,x, phi1, phi2,aHat,bHat] = gradientDescent(u,pole,gamma,a,b)
equationsSystem = @(t,y) prepareEquationsSystemToBeSolved(t,y,u,gamma,pole,a,b);
[tFromODE, y] = ode45(equationsSystem,[0 10],[0; 0; 0; 0; 0]);
x = y(:,1);
phi1 = y(:,2); phi2 = y(:,3);
aHat = -pole -y(:,4);
bHat = y(:,5);
end

function equations = prepareEquationsSystemToBeSolved(t,y,u,gamma,pole,a,b)
x = y(1); phi1 = y(2); phi2 = y(3); theta1 = y(4); theta2 = y(5);
xDot =  -a*x + b*u(t);
phi1Dot = pole*phi1 + x;
phi2Dot = pole*phi2 + u(t);
e = x - theta1*phi1 - theta2*phi2;
theta1Dot = gamma*e*phi1;
theta2Dot = gamma*e*phi2;

equations = [xDot; phi1Dot; phi2Dot; theta1Dot; theta2Dot];
end
