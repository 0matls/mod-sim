%%  System Modeling & Simulation
%   1st assignement, task 2, function 1 of 1

function theta = thetaEstimator2(pole1, pole2, samples, VcSamples, input1, input2, samplingTimes)

filter = [1 -(pole1+pole2) pole1*pole2];

phiMatrix = NaN(samples, 6); % phi-matrix with 6 parameters and #samples measurements
phiMatrix(:,1) = lsim(tf(-[1 0], filter), VcSamples, samplingTimes);
phiMatrix(:,2) = lsim(tf(-[0 1], filter), VcSamples, samplingTimes);
phiMatrix(:,3) = lsim(tf([1 0],filter), input1(samplingTimes), samplingTimes);
phiMatrix(:,4) = lsim(tf([0 1],filter), input1(samplingTimes), samplingTimes);
phiMatrix(:,5) = lsim(tf([1 0],filter), input2, samplingTimes);
phiMatrix(:,6) = lsim(tf([0 1],filter), input2, samplingTimes); 

theta = VcSamples * phiMatrix / (phiMatrix.' * phiMatrix);

end
