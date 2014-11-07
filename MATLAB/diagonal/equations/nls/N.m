function ns = N(~,yh,~,input_is_3d)
%L Returns nonliner function
if(input_is_3d)
    fft_dim = 3;
else
    fft_dim = 2;    
end
y = ifft(yh,[],fft_dim);
ns = 2*1i*fft(y.^2 .* conj(y),[],fft_dim);
end