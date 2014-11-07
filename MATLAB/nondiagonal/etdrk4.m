function [ts, ys, tcpu, tccpu] = etdrk4(L,N,tspan,y0,Nt,options)
% ETDRK4 Implements ETDRK4 Numerical Time Integrator
% PARAMETERS
%   L       - function, returns linear operator
%   N       - function, nonlinear operator
%   tspan   - integration bounds
%   y0      - initial condition
%   Nt      - number of timesteps
%   options - struct with fields: 
%               max_ts_to_store     - max solution values to store
%               problem_parameters  - struct which is passed to L and N Functions
% RETURNS
%   ys      - solution values
%   ts      - times cooresponging to solution values
%   tcpu    - seconds to compute solution
%   tccpu   - second to initialize coefficients

% === START Load Options ===
params  = options.parameters;
if(isfield(options,'max_ts_to_store'))
    max_ts_to_store = max(2,options.max_ts_to_store);
else
    max_ts_to_store = max(2,min(5000000/length(y0),1000));
end
% === END Load Options === 

% Numerical Parameters
h  = (tspan(end)-tspan(1))/Nt;
t0 = tspan(1);

% Data Parameters
skip_rate = ceil(Nt/ max_ts_to_store);
ys      = zeros(length(y0),ceil(Nt/skip_rate)+1);
ts      = zeros(size(ys,2),1);
ys(:,1) = y0; save_count = 2;
ts(1)   = tspan(1);

% === START Initialize Coefficients ====
tic;
[E, E2, A, b1, b23, b4] = initETDCoefficients(L,h);
tccpu = toc;
% === END Initialize Coefficients ====

tic;
for i=1:Nt  
    Nv = N(t0,y0,params);
    a  = E2*y0 + 0.5*A*Nv;
    Na = N(t0+h/2,a,params);
    b  = E2*y0 + 0.5*A*Na;
    Nb = N(t0+h/2,b,params);
    c  = E2*a + A*(Nb - 0.5*Nv);
    Nc = N(t0+h,c,params);
    y0 = E*y0 + b1*Nv + b23*(Na + Nb) + b4*Nc;
    t0 = t0 + h;
    % Save Data
    if(mod(i,skip_rate) == 0 || i == Nt)
        ys(:,save_count) = y0;
        ts(save_count,:) = t0;
        save_count = save_count + 1;
    end   
end
tcpu = toc;
end


function [E, E2, A, b1, b23, b4] = initETDCoefficients(Lambda,h)
%INITETDRK4 Initializes coefficients for ETDRK4 for matrix
    E = expm(h*Lambda);
    if(issparse(Lambda))
        E = sparse(E);
    end
    P = phi(h/2*Lambda,1);
        E2  = P{1};
        A   = h*P{2};
    P = phi(h*Lambda,4);
        b1  = h*(4*P{4} - 3*P{3} + P{2});
        b23 = h*(-4*P{4} + 2*P{3});
        b4  = h*(4*P{4} - P{3});
end