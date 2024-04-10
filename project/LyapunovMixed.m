clc; clear; close all;
%%  System Modeling & Simulation
% project - main script - 1 of 2 scripts 
%%
% a = ??;
% b = ??;
xDesired = @(t) 1*cos(3*t);  
gammaValues = [1 1];
thetaM = 2;
[tFromODE,x,xHat,aHat,bHat] = lyapunovMixed(gammaValues,thetaM,xDesired);
%%
rZero = 4;
rInf = 0.01; %fixed 
lambda = 2; %fixed
rho = @(t)  (rZero - rInf)*exp(-lambda*t) + rInf;
k = 2;
z = @(t,x) (x - xDesired(t)) ./ rho(t);
epsilonFunction = @(t,x) log( (1 + z(t,x)) / (1 - z(t,x)) ); 
r = @(t,x) 2 ./ rho(t) ./ (1 - z(t,x).^2 );
u = @(t,w,x) -k./w .* (r(t,x) + 1./r(t,x)).*epsilonFunction(t,x);
error = x - xHat;
error_x_x_d = x - xDesired(tFromODE);
%%
figure(1)
subplot(2,2,1)
plot(tFromODE, error_x_x_d)
xlabel('Time (s)','interpreter','latex')
title('$-\rho(t) < x(t) -x_{d}(t) < \rho(t) $','interpreter','latex');
hold on
plot(tFromODE, rho(tFromODE))
plot(tFromODE, -rho(tFromODE))
legend("x - x_{desired}", "rho(t)", "-rho(t)") 
subplot(2,2,2)
plot(tFromODE, x);
xlabel('Time (s)','interpreter','latex')
hold on 
plot(tFromODE, xHat);
plot(tFromODE, error);
legend("x", "xHat", "error") 
title('$x - \hat{x} \quad and \quad their \quad error $','interpreter','latex');
subplot(2,2,3)
plot(tFromODE, u(tFromODE,bHat,x))
xlabel('Time (s)','interpreter','latex')
title('$u(t) $','interpreter','latex');
subplot(2,2,4)
hold on
xlabel('Time (s)','interpreter','latex')
title('$\hat{a} \quad and \quad \hat{b}\quad values $','interpreter','latex')
plot(tFromODE, bHat);
plot(tFromODE, aHat);
legend("bHat", "aHat") 
%%
function [tFromODE,x,xHat,aHat,bHat] = lyapunovMixed(gamma1gamma2,thetaM,xDesired)
rZero = 4;
rInf = 0.01; %fixed 
lambda = 2; %fixed
rho = @(t)  (rZero - rInf)*exp(-lambda*t) + rInf;
k = 2;
z = @(t,x) (x - xDesired(t))/rho(t);
epsilonFunction = @(t,x) log( (1 + z(t,x)) / (1 - z(t,x)) ); 
r = @(t,x) 2 / rho(t) / (1 - z(t,x)^2 );
u = @(t,w,x) -k/w * (r(t,x) + 1/r(t,x))*epsilonFunction(t,x);
equationsSystem = @(t,y) prepareEquationsSystemToBeSolved_MixedLyapunov(t,y,gamma1gamma2,thetaM, u);
[tFromODE,y] = ode15s(equationsSystem, [0 50], [0; 0; 0.3; 0.2]);
x = y(:,1);
xHat = y(:,2);
aHat = y(:,3);
bHat = y(:,4);

end

function equations = prepareEquationsSystemToBeSolved_MixedLyapunov(t,y,gamma1gamma2,thetaM, u)
x = y(1); xHat=y(2); theta1 = y(3); theta2 = y(4);
w = theta2;

xDot = sys(x, u(t,w,x));
error = x - xHat;
xHatDot = theta1*x + theta2*u(t,w,x) + thetaM*error;
thetaDot1 = gamma1gamma2(1)*error*x;
thetaDot1 = projectionTheta1(theta1,thetaDot1, 0.001);

thetaDot2 = gamma1gamma2(2)*error*u(t,w,x);  
thetaDot2 = projectionTheta2(theta2, thetaDot2, 0.001);
equations = [xDot; xHatDot; thetaDot1; thetaDot2];
% t
end

function [thetaDot1] = projectionTheta1(theta1, thetaDot1, epsilon)
    if theta1>=(-epsilon) && theta1 < 0 && thetaDot1 < 0
        thetaDot1 = 1 - min(1,-theta1/epsilon)*thetaDot1;
    end
end

function [thetaDot2] = projectionTheta2(theta2, thetaDot2, epsilon)
  
   if (theta2 < 0.1)  && thetaDot2 < 0% && theta2 > (0.1 - 10*epsilon)
       thetaDot2 = (1 - min(1, (0.1-theta2)/epsilon) )*thetaDot2;  
   elseif theta2 > 2  && thetaDot2 > 0  && theta2 < (2+epsilon)
        thetaDot2 = (1 - min(1,(theta2 - 2)/epsilon))*thetaDot2;
    end
end