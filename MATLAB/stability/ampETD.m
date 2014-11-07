function psi = ampETD(n,m,r,z,I)
%AMP_ETD Computes amplification factor (labeled \psi(r,z) in paper) for ETD 
% spectral deferred correction scheme assuming n chebyshyshev quadrature 
% points and m corrections sweeps.
%PARAMETERS
% n  - number of quadrature points
% m  - number of correction sweeps
% r  - h*\lambda
% z  - h*\mu
% I  - ETD integration matrix

% Parameters
xj = 0.5 - 0.5*cos(pi*(0:n-1)/(n-1));
h = xj(2:end) - xj(1:end-1);

% Countour Integral Parameters
radius     = 2; 
num_points = 32;
min_r      = 1;

% Determine T1(i) = (exp(r*h(i)) - 1)/r
T1 = zeros(n-1,1);
for i=1:n-1
    if(abs(h(i)*r) < min_r)  % Contour Integral Approximation
        zr = h(i)*r + radius*exp(2*pi*1i*(0:(num_points-1))/num_points);
        T1(i) = h(i)*mean((exp(zr) - 1)./zr);
    else % Direct Formula    
        T1(i) = (exp(r*h(i)) - 1)/r;
    end
end

% Euler
psi = ones(n,1); 
for i=1:n-1
    psi(i+1) = (exp(r*h(i)) + z*T1(i))*psi(i);   
end

% Correction Sweeps
for k=1:m
    psi_old = psi;
    for i=1:n-1
        psi(i+1) = exp(r*h(i))*psi(i) + T1(i)*z*(psi(i) - psi_old(i)) + z*I(i,:)*psi_old;
    end
end
psi = psi(end);
end