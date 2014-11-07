function ns = N(~,yh,params,input_is_3d)
%L Returns nonliner function
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
ks = [0:Nx/2 -Nx/2+1:-1] * (2*pi/Lx);
if(input_is_3d)
    fft_dim = 3;
    ks = reshape(ks,[1 1 Nx]);
else
    fft_dim = 2;    
end
ns = -0.5i * bsxfun(@times,ks,fft(ifft(yh,[],fft_dim).^2,[],fft_dim));
end