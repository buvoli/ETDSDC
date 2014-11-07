function cs = L(params)
%L Returns coefficients for linear ODEs

%set up wavenumbers
wn_x = [0:params.Nx/2 -params.Nx/2+1:-1].' * (2*pi/(params.Lx));
wn_y = [0:params.Ny/2 -params.Ny/2+1:-1].' * (2*pi/(params.Ly));

% set up Linear L Matrix
m_x = spdiags(wn_x,0,params.Nx,params.Nx);
m_y = spdiags(wn_y,0,params.Nx,params.Nx);
O    = ones(params.Ny,params.Nx);

% Inverse Laplace Operator
LA  = ((1i*m_y).^2)*O + O*((1i*m_x).^2);
ILA = 1./LA;
ILA(1,1) = 0; %Set (0,0) mode to zero

% Linear Operator
L = (-params.beta*O*(1i*m_x)).*ILA - params.epsilon*O - params.v * ((1i*m_y).^8*O + O*(1i*m_x).^8);

% Reshape to row vector format
cs = reshape(L,1,params.Nx*params.Ny);
end