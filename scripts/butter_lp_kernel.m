function kl = butter_lp_kernel(N, Dl, no),
% butterworth low pass filter kernel 
% no: order , N : size of the matrix filter
D=false(N,N);
r=N/2;
D(round(r),round(r))=1;
D=bwdist(logical(D),'euclidean');
denom = 1+(sqrt(2)-1).*(Dl./D).^(2*no);
kl = 1-1./(denom);
%figure, imshow(kl,[]); title ('low pass butterw')
end