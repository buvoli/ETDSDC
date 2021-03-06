% Solves KDV equation in Fourier space (diagonal Lambda)

% Spatial Descretization
Lx = 2; Nx = 2^8;
xs = linspace(0,Lx,Nx+1)'; xs(end) = [];

% Temporal Discretization
tspan = [0 3.6/pi]; Nt = 500;

% initial conditions
y0 = cos(pi*xs);
pars = struct('Lx',Lx,'Nx',Nx,'delta',.022);

% set L and N
LF = L(pars);
NF = @N;

% Plot/Error Variables
error_filters = {@(x,y) max(abs(x-y))./max(abs(y))};
filter        = @(x) x;
equation_dimension = 1;