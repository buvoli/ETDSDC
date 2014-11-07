% Spatial Descretization
Lx = 2; Nx = 2^8;
xs = linspace(0,Lx,Nx+1)'; xs(end) = [];

% Temporal Discretization
tspan = [0 3.6/pi];

% initial conditions
y0 = fft(cos(pi*xs));
pars = struct('Lx',Lx,'Nx',Nx,'delta',.022);

% set L and N
LF = L(pars);
NF = @N;

% Plot/Error Variables
error_filter = @(x,y) max(abs(real(ifft(x-y))))./max(abs(real(ifft(y))));
filter       = @(x) real(ifft(x));
equation_dimension = 1;