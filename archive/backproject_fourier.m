function fourier_radial=backproject_fourier(f_p, prj_angles)
	interp_m='spline';

	% prj_angles = round(prj_angles*100)/100;

	[prj_angles, a_order] = sort(prj_angles);
	f_p = f_p(:, a_order);

	[prj_angles, idx] = unique(prj_angles, 'stable');
	f_p = f_p(:, idx);

	nfp=length(f_p(:,1)); % final length of the fft1 

	%% gridding for polar to cartesian 
	omega_sino=(-(nfp-1)/2:(nfp-1)/2).*(2*pi/size(f_p,1));
	theta=prj_angles*pi/180;
	omega_image=omega_sino;
	[theta_grid, omega_grid] = meshgrid(theta,omega_sino); 

	% grid image target 
	[omega_grid_x, omega_grid_y] = meshgrid(omega_image, omega_image);
	[coo_th_fft2, coo_r_fft2] = cart2pol(omega_grid_x,omega_grid_y);


	coo_r_fft2=coo_r_fft2.*sign(coo_th_fft2); % if theta >pi
	coo_th_fft2(coo_th_fft2<0)=coo_th_fft2(coo_th_fft2<0)+pi;

	fourier_radial = interp2(theta_grid,omega_grid,f_p,coo_th_fft2,coo_r_fft2,interp_m,0);
end