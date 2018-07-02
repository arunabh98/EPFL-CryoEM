function projection = project_fourier_alternate(fourier_radial, probe_angle, shift, proj_length)
	% Resolution of the grid.
	resolution_grid = 300;

	% The entire array of angles. 
	all_prj_angles = 0:0.1:179.9;

	% Length of the projections.
	nfp = proj_length;

	% Create the grid.
	omega_sino = (-(nfp-1)/2:(nfp-1)/2).*(2*pi/proj_length);
	all_theta = all_prj_angles*pi/180;
	[theta_grid, omega_grid] = meshgrid(all_theta, omega_sino);

	[x, y] = pol2cart(theta_grid, omega_grid);

	x = round(x*resolution_grid);
	y = round(y*resolution_grid);

	min_value = min(min(x(:)),  min(y(:)));

	probe_theta = probe_angle*pi/180;
	[probe_theta_grid, probe_omega_grid] = meshgrid(probe_theta, omega_sino);
	[probe_x, probe_y] = pol2cart(probe_theta_grid, probe_omega_grid);
	probe_x = round(probe_x*resolution_grid);
	probe_y = round(probe_y*resolution_grid);
	probe_x = probe_x - min_value + 1;
	probe_y = probe_y - min_value + 1;

	projection = zeros(proj_length, 1);

	parfor i=1:proj_length
		projection(i) = fourier_radial(probe_y(i), probe_x(i));
	end

	% Calculate the frequencies in the projection.
	first_part = (1/(nfp*2):1/nfp:0.5)';
	second_part = -first_part(1:end-1);
	frequencies = [first_part; flipud(second_part)];

	% Shift the projection.
    if mod(shift, 2) == 0
        freq_compensation = exp(-1j*2*pi*frequencies*shift);
    else
        freq_compensation = -exp(-1j*2*pi*frequencies*shift);
    end
    projection = projection.*freq_compensation;
end
