function fourier_radial=backproject_fourier_alternate(f_p, prj_angles)
	% Resolution of the grid.
	resolution_grid = 1000;

	% The entire array of angles. 
	all_prj_angles = 0:179;

	% Length of the projections.
	nfp=length(f_p(:,1));

	% Create the grid.
	omega_sino = (-(nfp-1)/2:(nfp-1)/2).*(2*pi/size(f_p,1));
	all_theta = all_prj_angles*pi/180;
	[theta_grid, omega_grid] = meshgrid(all_theta, omega_sino); 

	[x, y] = pol2cart(theta_grid, omega_grid);

	x = round(x*resolution_grid);
	y = round(y*resolution_grid);

	min_value = min(min(x(:)),  min(y(:)));
	x = x - min_value + 1;
	y = y - min_value + 1;

	size_matrix = max(max(x(:)), max(y(:)));
	fourier_radial = zeros(size_matrix, size_matrix);
	count_matrix = zeros(size_matrix, size_matrix);

	probe_theta = prj_angles*pi/180;
	[probe_theta_grid, probe_omega_grid] = meshgrid(probe_theta, omega_sino); 
	[probe_x, probe_y] = pol2cart(probe_theta_grid, probe_omega_grid);
	probe_x = round(probe_x*resolution_grid);
	probe_y = round(probe_y*resolution_grid);
	probe_x = probe_x - min_value + 1;
	probe_y = probe_y - min_value + 1;

	for j=1:size(prj_angles, 2)
		for i=1:nfp
			fourier_radial(probe_y(i, j), probe_x(i, j)) = fourier_radial(probe_y(i, j), probe_x(i, j)) + f_p(i, j);
			count_matrix(probe_y(i, j), probe_x(i, j)) = count_matrix(probe_y(i, j), probe_x(i, j)) + 1;
		end
	end

	index_non_nan = find(count_matrix ~= 0);
	fourier_radial(index_non_nan) = ...
		fourier_radial(index_non_nan)./count_matrix(index_non_nan);

	fourier_radial(fourier_radial == 0) = NaN;
	fourier_radial = inpaint_nans(fourier_radial, 2);
end