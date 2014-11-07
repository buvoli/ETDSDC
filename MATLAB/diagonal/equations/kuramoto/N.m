function ns = N(~,yh,params,input_is_3d)
%L Returns nonliner function
Lx = params.Lx;     % Lx - domain size
Nx = params.Nx;     % Nx - number of spatial points
if(input_is_3d)
    ks = reshape([0:Nx/2 -Nx/2+1:-1] * (2*pi/Lx),1,1,Nx);
    ns = -0.5i * bsxfun(@times,ks,fft(ifft(yh,[],3).^2,[],3));
else
    ks = [0:Nx/2 -Nx/2+1:-1] * (2*pi/Lx);
    ns = -0.5i * bsxfun(@times,ks,fft(ifft(yh,[],2).^2,[],2));
end
end