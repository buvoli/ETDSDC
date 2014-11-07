function [W, P01]=initW(L,h,tau)
%INITW_SCALER Initializes ETDSDC W functions for scaler/vector L using weights function by Fornberg.
% Input Parameters
%   L   - (array)  L = [\Lambda_1, \Lambda_2, ..., \Lambda_N]
%   h   - (double) timestep
%   tau - (array)  normalized quadrature points
% Output Parameters
%   W   - (array) array of dimensions (n-1) x n x length(L). 
%         W(:,:,j) contains integration matrix cooresponding to L(j)
%   P01 - (array) array containg phi_0 and phi_1 needed for ETD Euler method.

eta = tau(2:end) - tau(1:end-1);
n   = length(tau);
W   = zeros(n-1,n,length(L));             % stores integration matrices W
P01 = zeros(2,length(L),n-1);             % stores phi_0 and phi_1 for ETD Euler
for i=1:n-1
    q = (tau - tau(i))/eta(i);            % scaled quadrature ppints
    a = weights(0,q,n-1);                 % finite difference matrix
    p = phi(eta(i)*h*L,n);         % phi functions 0 to n          
    W(i,:,:) = h*eta(i)*a*p(2:end,:);     % store ith row of W matrix
    P01(1,:,i) = p(1,:);                  % store phi_0 and phi_1
    P01(2,:,i) = h*eta(i)*p(2,:);         % store phi_0 and phi_1
end
end