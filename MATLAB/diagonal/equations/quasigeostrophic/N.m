function ns = N(~,yh,params,input_is_3d)
%N Returns nonliner function
% yh - correponds to \del^2 \phi

num_yh = size(yh,1);
if(input_is_3d)
    del2phi_h = reshape(permute(yh,[3 2 1]),params.Ny,params.Nx,num_yh);
else
    del2phi_h = reshape(permute(yh,[2 1]),[params.Ny,params.Nx,num_yh]);
end

del2phi = ifft2(del2phi_h);

DX = params.DX;
DY = params.DY;
IL = params.IL; 

phi_h = bsxfun(@times,IL,del2phi_h);
N1 = bsxfun(@times,DX,fft2(bsxfun(@times,ifft2(bsxfun(@times,-DY,phi_h)),del2phi)));
N2 = bsxfun(@times,DY,fft2(bsxfun(@times,ifft2(bsxfun(@times,DX,phi_h)),del2phi)));
N  = -(N1 + N2);

if(input_is_3d)
    ns = permute(reshape(N,[params.Nx*params.Ny 1 num_yh]),[3 2 1]);
else
    ns = permute(reshape(N,[params.Nx*params.Ny,num_yh]),[2 1]);
end
end

% NOTE:

% EQUIVALENT LOOP CODE (for converting 3D vector into appropriate form
% for fft2)
% del2phi_h = zeros(params.Ny,params.Nx,num_yh);
% for i=1:num_yh
%     del2phi_h(:,:,i) = reshape(yh(i,1,:),params.Ny,params.Nx); 
% end

% EQUIVALENT LOOP CODE (for convering vector back to 3D
%ns = zeros(num_yh,1,params.Ny*params.Nx);
%for i=1:num_yh
%    ns(i,1,:) = reshape(N(:,:,i),1,params.Nx*params.Ny);
%end