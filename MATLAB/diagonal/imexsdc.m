function [ts,ys,tcpu,tcpuc] = imexsdc(L,N,tspan,y0,Nt,options)
% IMEXSDC Implements IMEX Spectral Deferred Correction Scheme
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
%               quadpts             - nx1 array which contains normalized
%                                     quadrature points
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
h   = (tspan(end)-tspan(1))/Nt;
eta = tau(2:end) - tau(1:end-1);
y_stage = zeros(n,length(y0));
yn  = zeros(n,length(y0));

% Data Parameters
skip_rate = ceil(Nt/ max_ts_to_store);
ys   = zeros(length(y0),ceil(Nt/skip_rate)+1);
ts   = zeros(size(ys,2));
ys(:,1) = y0; save_count = 2;
ts(1)   = tspan(1); t0 = tspan(1);

% === START Initialize Coefficients ====
tic;
L  = reshape(L,1,length(L));
D  = (1 - kron(h*eta,L)); 
IM = initIM(h,tau);
tcpuc = toc;
% === END Initialize Coefficients ====

tic;
y_stage(1,:) = y0; % store solution at each stage in nxm matrix where A(n,m) is solution at \tau_n, x_m
for i=1:Nt
    %sub stage
    for j=1:n-1
        yn(j,:) = N(t0+h*tau(j),y_stage(j,:),params,false);
        y_stage(j+1,:) = (y_stage(j,:) + h*eta(j)*yn(j,:))./D(j,:);
    end
    yn(n,:) = N(t0+h*tau(n),y_stage(n,:),params,false);
    for k=1:m
        I = IM*(repmat(L,n,1).*y_stage + yn);
        y_stage(2,:) = (y_stage(1,:) - h*eta(1)*L.*y_stage(2,:) + I(1,:))./D(1,:);
        for j=2:n-1
            n_new = N(t0 + h*tau(j),y_stage(j,:),params,false);
            y_stage(j+1,:) = (y_stage(j,:) - h*eta(j)*L.*y_stage(j+1,:) + h*eta(j)*(n_new - yn(j,:)) + I(j,:))./D(j,:);
            yn(j,:) = n_new;
        end
        yn(n,:) = N(t0+h*tau(n),y_stage(n,:),params,false);
    end
    y_stage(1,:) = y_stage(end,:);
    t0 = t0 + h;
    % Save Data
    if(mod(i,skip_rate) == 0 || i == Nt)
        ys(:,save_count) = y_stage(1,:);
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