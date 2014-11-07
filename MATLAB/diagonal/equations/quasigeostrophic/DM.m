function [DX, DY, IL, DEL2] = DM(params)
%DM Returns differentation Matrices & Inverse Lapalacian for 2D periodic Domain

%set up wavenumbers
wn_x = [0:params.Nx/2 -params.Nx/2+1:-1].' * (2*pi/(params.Lx));
wn_y = [0:params.Ny/2 -params.Ny/2+1:-1].' * (2*pi/(params.Ly));

% set up Linear L Matrix
m_x = spdiags(wn_x,0,params.Nx,params.Nx);
m_y = spdiags(wn_y,0,params.Nx,params.Nx);
O   = ones(params.Ny,params.Nx);

% Build Differentiation coefficients
DX  = O*(1i*m_x);
DY  = (1i*m_y)*O;

% Build Invese Laplace coefficients 
DEL2    = (1i*m_y).^2*O + O*(1i*m_x).^2;
IL      = 1./DEL2;
IL(1,1) = 0; %Set (0,0) mode to zero

end