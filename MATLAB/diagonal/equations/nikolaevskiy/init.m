% Spatial Descretization
Lx = 150*pi; Nx = 2^12;
xs = linspace(-Lx/2,Lx/2,Nx+1); xs(end) = [];

% Temporal Discretization
tspan = [0 50];

% initial conditions
y0 = fft(sin(xs) + 0.1*sin(xs/25)); 
pars = struct('Lx',Lx,'Nx',Nx,'r',1/4,'alpha',2.1,'beta',0.77);

% set L and N
LF = L(pars);
NF = @N;

% Plot/Error Variables
error_filter = @(x,y) max(abs(real(ifft(x-y))))./max(abs(real(ifft(y))));
filter       = @(x) real(ifft(x));
equation_dimension = 1;