function ns = N(~,y,params)
%L Returns nonliner function
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1]' * (2*pi/Lx); ks(Nx/2+1) = 0;
ns = (ifft(repmat(-1i*ks,1,size(y,2)).*fft(y.^2)/2));
end