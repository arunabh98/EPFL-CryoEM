function fourier_radial=backproject_fourier_alternate(f_p, prj_angles, shifts)
	% Length of the projections.
	nfp=length(f_p(:,1));

	% Calculate the frequencies in the projection.
	first_part = (1/(nfp*2):1/nfp:0.5)';
	second_part = -first_part(1:end-1);
	frequencies = [first_part; flipud(second_part)];

	% Shift back the projections.
	parfor i=1:size(f_p, 2)
	    shift = shifts(i);
	    if mod(shift, 2) == 0
	        freq_compensation = exp(1j*2*pi*frequencies*shift);
	    else
	        freq_compensation = -exp(1j*2*pi*frequencies*shift);
	    end
	    f_p(:, i) = f_p(:, i).*freq_compensation;
	end

	% Resolution of the grid.
	resolution_grid = 300;

	% The entire array of angles. 
	all_prj_angles = 0:0.1:179.9;

	% Create the grid.
	omega_sino = (-(nfp-1)/2:(nfp-1)/2).*(2*pi/size(f_p, 1));
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

	[prj_angles, prj_sort] = sort(prj_angles);
	f_p = f_p(:, prj_sort);
	[unique_angles, unique_indices] = unique(prj_angles);

	for i=2:size(unique_indices, 1) - 1
		same_projections = f_p(:, unique_indices(i):unique_indices(i+1) - 1);
		prev_projection_distance_matrix =...
			bsxfun(@minus, same_projections, f_p(:, unique_indices(i) - 1));
		[~, order_proj] = sort(vecnorm(prev_projection_distance_matrix));
		f_p(:, unique_indices(i):unique_indices(i+1) - 1) =...
            same_projections(:, order_proj);
        num_projections = size(same_projections, 2);
        left_half = unique_angles(i):0.5:(unique_angles(i) +...
            0.5*floor(num_projections/2));
        right_half = (unique_angles(i) - 0.5*round(num_projections/2) +...
           0.5):0.5:(unique_angles(i) - 0.5);
        corr_angles = [right_half left_half];
        prj_angles(:, unique_indices(i):unique_indices(i+1) - 1) = corr_angles;
	end

	probe_theta = prj_angles*pi/180;
	[probe_theta_grid, probe_omega_grid] = meshgrid(probe_theta, omega_sino); 
	[probe_x, probe_y] = pol2cart(probe_theta_grid, probe_omega_grid);
	probe_x = round(probe_x*resolution_grid);
	probe_y = round(probe_y*resolution_grid);
	probe_x = probe_x - min_value + 1;
	probe_y = probe_y - min_value + 1;

	for j=1:size(prj_angles, 2)
		for i=1:nfp
			fourier_radial(probe_y(i, j), probe_x(i, j)) = ...
                fourier_radial(probe_y(i, j), probe_x(i, j)) + f_p(i, j);
			count_matrix(probe_y(i, j), probe_x(i, j)) = ...
                count_matrix(probe_y(i, j), probe_x(i, j)) + 1;
		end
	end

	index_non_nan = find(count_matrix ~= 0);
	fourier_radial(index_non_nan) = ...
		fourier_radial(index_non_nan)./count_matrix(index_non_nan);

	fourier_radial(fourier_radial == 0) = NaN;
	fourier_radial = inpaint_nans(fourier_radial, 5);
end
