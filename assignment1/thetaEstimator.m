%%  System Modeling & Simulation
%   1st assignement, task 1, function 1 of 2

function theta = thetaEstimator(pole1, pole2, y, u, timespan)

filter = [1 -(pole1+pole2) pole1*pole2]; % to create Lambda(s) coefficients

phiMatrix = NaN(length(timespan), 3); % phi-matrix with 3 parameters and #timespan measurements
phiMatrix(:,1) = lsim(tf(-[1 0], filter), y(:,1), timespan);
phiMatrix(:,2) = lsim(tf(-1, filter), y(:,1), timespan);
phiMatrix(:,3) = lsim(tf(1, filter), u(timespan), timespan);

if rank(phiMatrix' * phiMatrix) < 3
    fprintf("The matrix is not inversible");
else
theta = y(:,1)' * phiMatrix / (phiMatrix' * phiMatrix);
end

end
