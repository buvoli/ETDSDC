function IM=initIM(h,tau)
%INITW_SCALER Returns (n-1)x(n) Integration matrix that integrates from tau(i) to tau(i+1)
% Parameters
%   tau - (vector) normalized quadrature points
% Output
%   I   - integration matrix
n   = length(tau);
eta = tau(2:end) - tau(1:end-1);
IM  = zeros(n-1,n);
p   = (1./factorial(1:n)).';          % phi functions for lambda=0
for i=1:n-1
	q  = (tau - tau(i))/eta(i);     % scaled quadrature ppints
	FD = weights(0,q,n-1);          % finite difference matrix
    IM(i,:) = h*eta(i)*FD*p;        % form ith row
end
end