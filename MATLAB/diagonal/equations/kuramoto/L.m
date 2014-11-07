function cs = L(params)
%L Returns coefficients for linear ODEs
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1] * (2*pi/Lx);
cs = -1 * params.epsilon * ((1i*ks).^2 + (1i*ks).^4);% - params.delta * (1i*ks).^3;
end

