function I_a = Ifft2_2_Img(I_fft2, crop)
    sz=size(I_fft2,1); mask=ones(sz,sz);
    Dl=0.6*sz;
    kl = butter_lp_kernel(sz, Dl, 2);
    I_fft2=I_fft2.*kl;
    target=fftshift(ifft2(ifftshift(I_fft2))); 

    target=(target-min(target(:)))./(max(target(:))-min(target(:)));

    if crop>0
      
        target=target(crop:end-crop,crop:end-crop);
        
    end
    I_a= abs(target); % module 
end