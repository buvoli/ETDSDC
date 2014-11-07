function [ts,ys,tcpu,tcpuc] = etdsdc(L,N,tspan,y0,Nt,options)
% ETDSDC Implements ETDSDC Scheme
% PARAMETERS
%   L       - (array) linear operator
%   N       - function, nonlinear operator
%   tspan   - integration bounds
%   yi      - initial condition
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
if(isfield(options,'quadpts'))
    tau = options.quadpts;
else
    tau  = 0.5 - 0.5*cos(pi*((0:n-1)')/(n-1));
end
% === END Load Options === 

% Numerical Parameters
h    = (tspan(end)-tspan(1))/Nt; 
t0 = tspan(1);
phi  = zeros(length(y0),n);
phi_n = zeros(length(y0),n);
I     = zeros(length(y0),n);

% Data Parameters
skip_rate = ceil(Nt/ max_ts_to_store);
ys   = zeros(length(y0),ceil(Nt/skip_rate)+1);
ts   = zeros(size(ys,2),1);
ys(:,1) = y0; save_count = 2;
ts(1)   = tspan(1);

% === START Initialize Coefficients ====
tic;
[W, P01] = initW(L,h,tau);
tcpuc = toc;
% === END Initialize Coefficients ====

tic;
phi(:,1) = y0;
for j=1:Nt
    % ETD Euler
    for i=1:n-1
        phi_n(:,i)  = N(t0 + h*tau(i),phi(:,i),params);
        phi(:,i+1)  = P01{1,i}*phi(:,i) + P01{2,i}*phi_n(:,i);
    end
    phi_n(:,n) = N(t0 + h*tau(n),phi(:,n),params);
    % Correction Sweeps
    for k=1:m
        % Compute Picard Iteration Terms (could eventually be done in parallel)
        for i=1:n-1
            I(:,i) = W{1,i}*phi_n(:,1);
            for l=2:n
                I(:,i) = I(:,i) + W{l,i}*phi_n(:,l);
            end
        end
        % Correction i = 1
        phi(:,2) = P01{1,1}*phi(:,1) + I(:,1);
        % Correction i = 2...n-1
        for i=2:n-1
            new_n = N(t0 + h*tau(i),phi(:,i),params);
            phi(:,i+1) = P01{1,i}*phi(:,i) + P01{2,i}*(new_n - phi_n(:,i)) + I(:,i);
            phi_n(:,i) = new_n;
        end
        phi_n(:,n) = N(t0 + h*tau(n),phi(:,n),params);
    end
    phi(:,1) = phi(:,n);
    t0 = t0 + h;
    % Save Data
    if(mod(j,skip_rate) == 0 || j == Nt)
        ys(:,save_count) = phi(:,1);
        ts(save_count,:) = t0;
        save_count = save_count + 1;
    end
end
tcpu = toc;
end


function [W, P01]  = initW(L,h,tau)
%INITW_MATRIX Initializes ETDSDC W functions for matrix A using weights function by Fornberg.
% PARAMETERS
%   L   - (matrix)  matrix corresponding to PDE linear operator
%   h   - (double)  timestep
%   tau - (array)   normalized quadrature nodes
% RETURNS
%   W   - (cell array) W{n,i} contains w_{n,i}(h*eta_i*A) for n=1,N
%   P   - (cell array) P{n,i} contains (h*eta_i)*Phi_{n,i}(h*eta_i*A) for n=0,1 
%                      which are needed for euler step.

% define the chebyshev points on I=[0 1]
n = length(tau);
eta = tau(2:n) - tau(1:n-1);
% Compute W Functions
P01 = cell(2,n-1);
W = cell(n,n-1);

for i=1:n-1
    P = phi(L*h*eta(i),n);   % init phi functions
    q = (tau - tau(i))/eta(i);                  % scaled quadrature ppints
    a = weights(0,q,n-1);                       % finite difference weights    
    for j=1:n
        W{j,i} = a(j,1)*P{2};
        for k=2:n
            W{j,i} = W{j,i} + a(j,k)*P{k+1}; % sum taylor series for each Tau Function
        end
        W{j,i} = h*eta(i)*W{j,i};
    end
    % Store P01
    P01{1,i} = P{1};
    P01{2,i} = h*eta(i)*P{2};
end
end