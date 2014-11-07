function psi = ampIMEX(n,m,r,z,I)
%AMP_IMEX Computes amplification factor (labeled \psi(r,z) in paper) for IMEX
% spectral deferred correction scheme assuming n chebyshyshev quadrature 
% points and m corrections sweeps.
%PARAMETERS
% n  - number of quadrature points
% m  - number of correction sweeps
% r  - h*\lambda
% z  - h*\mu
% I  - integration matrix

% Parameters
xj = 0.5 - 0.5*cos(pi*(0:n-1)/(n-1));
h = xj(2:end) - xj(1:end-1); 

% Euler's Method
psi = ones(n,1);
for i=1:n-1
    psi(i+1) = psi(i)*(1+z*h(i))/(1-r*h(i));
end

% Correction Sweeps
for k=1:m
    phi_old = psi;
    for i=1:n-1
        psi(i+1) = (psi(i) + h(i)*z*(psi(i) - phi_old(i)) - r*h(i)*phi_old(i+1) + I(i,:)*((r+z)*phi_old))/(1 - r*h(i));
    end
end
psi = psi(end);
end