%%  System Modeling & Simulation
%   1st assignement, task 1, function 2 of 2

function ds = msdODE(s, t, k, b, m, u) 
ds(1) = s(2);
ds(2) = -(b/m)*s(2) -(k/m)*s(1) + u(t)/m; 
ds = ds';
end