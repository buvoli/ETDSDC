function ns = N(~,yh,params)
%L Returns nonliner function
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1]' * (2*pi/Lx);
ns = repmat(-1i*ks,1,size(yh,2)).*(fft(ifft(yh).^2))/2;
end