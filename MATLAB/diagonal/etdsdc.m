function [ts,ys,tcpu,tcpuc] = etdsdc(L,N,tspan,y0,Nt,options)
% ETDSDC Implements ETDSDC Scheme
% PARAMETERS
%   L       - function, returns linear operator
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
h    = (tspan(end)-tspan(1))/Nt;
Nx   = length(y0);
y_stage = zeros(n,1,Nx);
yn = zeros(n,1,Nx);

% Data Parameters
skip_rate = ceil(Nt/ max_ts_to_store);
ys   = zeros(length(y0),ceil(Nt/skip_rate)+1);
ts   = zeros(size(ys,2),1);
ys(:,1) = y0; save_count = 2;
ts(1)   = tspan(1); t0 = tspan(1);

% === START Initialize Coefficients ====
tic;
[IM,P01] = initW(L,h,tau);
    % permute matrix IM and P01 so that we can use bsxfun
    IM  = permute(IM,[2 1 3]); 
    P01 = permute(P01,[3 1 2]); 
tcpuc = toc;
% === END Initialize Coefficients ====

tic;
y_stage(1,1,:) = y0; % store solution at each stage  in nx1xm matrix where A(i,1,m) is solution at \tau_n, x_m
for i=1:Nt
    %sub stage
    for j=1:n-1
        yn(j,1,:) = N(t0 + h*tau(j),y_stage(j,1,:),params,true);
        y_stage(j+1,1,:) = P01(j,1,:).*y_stage(j,1,:) + P01(j,2,:).*yn(j,1,:);
    end
    yn(n,1,:) = N(t0 + h*tau(n),y_stage(n,1,:),params,true);
    for k=1:m
        I = reshape(sum(bsxfun(@times,yn,IM)),n-1,1,Nx); % faster way of computing matrix product IM(:,:,i)*yn(:,:,i) for i=1...m
        y_stage(2,1,:) = P01(1,1,:).*y_stage(1,1,:) + I(1,1,:);
        for j=2:n-1
            n_new = N(t0+h*tau(j),y_stage(j,1,:),params,true);
            y_stage(j+1,1,:) = P01(j,1,:).*y_stage(j,1,:) + P01(j,2,:).*(n_new - yn(j,1,:)) + I(j,1,:);
            yn(j,1,:) = n_new;
        end
        yn(n,1,:) = N(t0+h*tau(n),y_stage(n,1,:),params,true);
    end
    y_stage(1,1,:) = y_stage(end,1,:);
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