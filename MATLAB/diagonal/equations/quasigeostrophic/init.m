% Spatial Descretization
Lx = 2*pi; Nx = 2^8;
Ly = 2*pi; Ny = 2^8;
xs = linspace(-Lx/2,Lx/2,Nx+1); xs(end) = [];
ys = linspace(-Ly/2,Ly/2,Nx+1); ys(end) = [];
[X, Y] = meshgrid(xs,ys);

% Temporal Discretization
tspan = [0 5];

% Parameters & Differentiation Matrices
pars = struct('Lx',Lx,'Nx',Nx,'Ly',Ly,'Ny',Ny,'epsilon',.01,'v',1e-14,'beta',10);
[DX, DY, IL, DEL2] = DM(pars); pars.DX = DX; pars.DY = DY; pars.IL = IL;

% initial conditions
y0 = DEL2.*(fft2((1/8) * exp(-8*(2*Y.^2 + 0.5*X.^2 - pi/4).^2)));
    y0 = reshape(y0,1,Nx*Ny);

% set L and N
LF = L(pars);
NF = @N;

% Plot/Error Variables
to_physical   = @(x) real(reshape(ifft2(reshape(x,pars.Ny,pars.Nx)),1,pars.Ny*pars.Nx));
error_filter = @(x,y) max(abs(to_physical(x-y)))./max(abs(to_physical(y)));
filter        = @(x) real(ifft2(reshape(x,pars.Ny,pars.Nx)));
equation_dimension = 2;