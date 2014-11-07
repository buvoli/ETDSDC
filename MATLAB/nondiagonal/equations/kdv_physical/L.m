function LM = L(params)
%L Returns coefficients for linear ODEs
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1]' * (2*pi/Lx); ks(Nx/2 + 1) = 0; % zero out mode N/2+1 (to insure odd derivative is real)
LM = real(ifft(eye(Nx)) * diag(-1*(params.delta)^2*(1i*ks).^3) * fft(eye(Nx)));
end