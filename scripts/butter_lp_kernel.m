function kl = butter_lp_kernel(N, Dl, no),
	D=false(N,N);
	r=N/2;
	D(round(r),round(r))=1;
	D=bwdist(logical(D),'euclidean');
	denom = 1+(sqrt(2)-1).*(Dl./D).^(2*no);
	kl = 1-1./(denom);
end