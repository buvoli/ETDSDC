function [ts,ys,tcpu,tcpuc] = imexsdc(L,N,tspan,y0,Nt,options)
% IMEXSDC Implements IMEXSDC scheme
% PARAMETERS
%   L       - vector, contains diagonal components of linear operator
%   N       - function, nonlinear operator
%   tspan   - integration bounds
%   y0      - initial condition
%   Nt      - number of timesteps
%   options - struct with fields: 
%               n                   - order of method
%               m                   - number of correction sweeps
%               max_ts_to_store     - max solution values to store
%               problem_parameters  - struct which is passed to L and N Functions
% RETURNS
%   ys      - solution values
%   ts      - times cooresponging to solution values
%   tcpu    - seconds to compute solution
%   tccpu   - second to initialize coefficients

% === START Load Options ===
n       = options.n;
m       = options.m;
params  = options.parameters;
if(isfield(options,'max_ts_to_store'))
    max_ts_to_store = max(2,options.max_ts_to_store);
else
    max_ts_to_store = max(2,min(5000000/length(y0),1000));
end
% === END Load Options === 

% Numerical Parameters
h     = (tspan(end)-tspan(1))/Nt; 
t0    = tspan(1);
tau   = 0.5 - 0.5*cos(pi*((0:n-1)')/(n-1));
phi   = zeros(length(y0),n);
phi_n = zeros(length(y0),n);

% Data Parameters
skip_rate = ceil(Nt/ max_ts_to_store);
ys        = zeros(length(y0),ceil(Nt/skip_rate)+1);
ts        = zeros(size(ys,2),1);
ys(:,1)  = y0; save_count = 2;
ts(1)     = tspan(1);

% === START Initialize Coefficients ====
tic;
% Lu factorization
eta = tau(2:end) - tau(1:end-1); 
Ls  = cell(n); Us = cell(n);
for i=1:n-1
    [l, u] = lu(speye(length(y0)) - h*eta(i)*L);
    Ls{i} = l; Us{i} = u;
end
% initialize integration matrix (assuming chebyshev points)
IM = transpose(initIM(h,tau));
tcpuc = toc;
% === END Initialize Coefficients ====

tic;
phi(:,1) = y0;
for i=1:Nt
    %sub stage
    for j=1:n-1
        phi_n(:,j) = N(t0+h*tau(j),phi(:,j),params);
        phi(:,j+1) = Us{j}\(Ls{j}\(phi(:,j) + h*eta(j)*phi_n(:,j)));
    end
    phi_n(:,n) = N(t0+h*tau(n),phi(:,n),params);
    for k=1:m
        I = (L*phi + phi_n)*IM;
        phi(:,2) = Us{1}\(Ls{1}\(phi(:,1) - h*eta(1)*L*phi(:,2) + I(:,1)));
        for j=2:n-1
            n_new = N(t0 + h*tau(j),phi(:,j),params);
            phi(:,j+1) = Us{j}\(Ls{j}\(phi(:,j) - h*eta(j)*L*phi(:,j+1) + h*eta(j)*(n_new - phi_n(:,j)) + I(:,j)));
            phi_n(:,j) = n_new;
        end
        phi_n(:,n) = N(t0+h*tau(n),phi(:,n),params);
    end
    phi(:,1) = phi(:,n);
    t0 = t0 + h;
    % Save Data
    if(mod(i,skip_rate) == 0 || i == Nt)
        ys(:,save_count) = phi(:,1);
        ts(save_count,:) = t0;
        save_count = save_count + 1;
    end
end
tcpu = toc;
end

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

% Call using
%IM = transpose((h/2) * initCI(n)); 
%
% function [M, DCTP]=initCI(n)
% %BUILDIM Returns n-1xn Integration matrix that integrates from x(i) to x(i+1)
% % for i=1,...,n-1 where x(i) is the ith chebyshev points. For derivation see
% % Norris, Gordon F. "Spectral integration and the numerical solution of
% % two-point boundary value problems." (1999). p.8-12
% %PARAMETERS:
% %   n - number of chebyshev points
% 
% % Init Permutation Matrix (so that ordering is -1 ... 1)
% PM = fliplr(eye(n)); 
% 
% % Init DCT Matrix
% DCT = zeros(n,n);
% c = [2 ones(1,n-2) 2];
% for i=1:n
%     for j=1:n
%         DCT(i,j) = (2/((n-1)*c(j)*c(i)))*cos((i-1)*(j-1)*pi/(n-1));
%     end
% end
% 
% % Init Integration Matrix
% if(n>2)
%     ud = [1/4 -1./(2*(1:n-3)) -1/((n-1)^2-1)];
%     ld = [1   1./(2*(2:n-1))];
% else
%    ud = 1/2;
%    ld = 1;
% end
% A = diag(ld,-1) + diag(ud,1);
% 
% % Init IDCT Matrix
% IDCT = zeros(n,n);
% for i=1:n
%     for j=1:n
%         IDCT(i,j) = cos((i-1)*(j-1)*pi/(n-1));
%     end
% end
% 
% % Projector
% P = zeros(n-1,n);
% for i=1:n-1
%     P(i,i) = 1;
%     P(i,i+1) = -1;
% end
% 
% % Form M Matrix
% M = fliplr(eye(n-1))*P*IDCT*A*DCT*PM;
% DCTP = DCT*PM;
% 
% end