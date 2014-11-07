% Spatial Descretization
Lx = 32*pi; Nx = 2^10;
xs = linspace(0,Lx,Nx+1)'; xs(end) = [];

% Temporal Discretization
tspan = [0 60];

% initial conditions
y0 = fft(cos(xs/16).*(1 + sin(xs/16)));
pars = struct('Lx',Lx,'Nx',Nx,'epsilon',1,'delta',.05);

% set L and N
LF = L(pars);
NF = @N;

% Plot/Error Variables
error_filter = @(x,y) max(abs(real(ifft(x-y))))./max(abs(real(ifft(y))));
filter       = @(x) real(ifft(x));
equation_dimension = 1;