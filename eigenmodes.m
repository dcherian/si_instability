% I have an ODE BVP problem for which I need to find the growth rate
% idea is to iterate over Ri, k and find eigenvector (mode) and eigenvalue (growth rate)
% using MATLAB's ODE solvers

% loop stuff

% write ODE as 2 equation system



%%

syms s

Ri=0.4;
a = -5*pi/180;
k = 1e-4;

z = linspace(0,-1,60);

AA = zeros(length(z),length(z));

for i=1:length(z)
    p = s - k*z(i)*sin(a);
    
    
end

