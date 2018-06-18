function I_a = Ifft2_2_Img(I_fft2, crop)

if (1) %  (1) cuts high frequencies & mitigates noise in the reconstruction
    sz=size(I_fft2,1); mask=ones(sz,sz);
    % one of the following filterng strategies 
    %%
    if 1,
    Dl=0.6*sz;
    kl = butter_lp_kernel(sz, Dl, 2);
    I_fft2=I_fft2.*kl;
    end
    %%
    
end
target=fftshift(ifft2(ifftshift(I_fft2))); 
%end

target=(target-min(target(:)))./(max(target(:))-min(target(:))); % nolrmalize [0 1]
%target = imadjust(target,[min(target(:)); max(target(:))],[0; 1])

if crop>0
  
    target=target(crop:end-crop,crop:end-crop);
    
end
I_a= abs(target); % module 
end