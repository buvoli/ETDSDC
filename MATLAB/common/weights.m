function C = weights(z,x,m)
% WEIGHTS Compute coefficients for finite difference approximation for the
% derivatives 1 to m at point z assuming data is known at points in array x.
% Based on the program "weights" in 
%   B. Fornberg, "Calculation of weights in finite difference formulas",
%   SIAM Review 40 (1998), pp. 685-691.
%
% Input Parameters
%	z - (double) location where approximations are to be accurate,
%	x - (array)  interpolation points
%	m - (int)    highest derivative for which weights are sought
% Output Parameters
%	C - (matrix) weights at grid locations x for derivative of order j<=m are found in c(:,j)

c1 = 1;
c4 = x(1) - z;
C = zeros(length(x),m+1); % coefficient array
C(1,1) = 1;

n = length(x);
for i=2:n
    mn = min(i,m+1);
    c2 = 1;
    c5 = c4;
    c4 = x(i) - z;
    for j=1:i-1
        c3 = x(i) - x(j);
        c2 = c2*c3;
        if(j == i-1)
            for k=mn:-1:2
                C(i,k) = c1*((k-1)*C(i-1,k-1) - c5*C(i-1,k))/c2;
            end
            C(i,1) = -c1*c5*C(i-1,1)/c2;
        end
        for k=mn:-1:2
            C(j,k) = (c4*C(j,k) - (k-1)*C(j,k-1))/c3;
        end
        C(j,1) = c4*C(j,1)/c3;
    end
    c1 = c2;
end
end