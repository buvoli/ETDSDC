% Spatial Descretization
Lx = 8*pi; Nx = 2^10;
xs = linspace(-Lx/2,Lx/2,Nx+1); xs(end) = [];

% Temporal Discretization
tspan = [0 15];

% initial conditions
y0 = fft(1 + (1/100)*exp(2*pi*1i*xs/Lx));
pars = struct('Lx',Lx,'Nx',Nx);

% set L and N
LF = L(pars);
NF = @N;

% Solution filter
error_filter = @(x,y) max(abs(ifft(x-y)))./max(abs(ifft(y)));
filter       = @(x) abs(ifft(x));
equation_dimension = 1;