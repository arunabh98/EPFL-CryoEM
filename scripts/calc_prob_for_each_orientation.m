function [prob_matrix, dist_matrix] = calc_prob_for_each_orientation(...
    f_image_estimate, f_projections, theta_estimate, first_estimate_theta,...
    first_estimate_shifts, shift_estimate,...
    noise_estimate, projection_parameters, prior_parameters)

	% Define prior parameters.
    max_angle_err = prior_parameters.max_angle_err;
    max_shift_err = prior_parameters.max_shift_err;
    resolution_angle = prior_parameters.resolution_angle;
    resolution_space = prior_parameters.resolution_space;

    % Calculate the probabilities of each orientation for each projection.
	prob_matrix_height = (2*max_angle_err)/resolution_angle + 1;
	prob_matrix_width = 2*max_shift_err/resolution_space + 1;
	prob_matrix = ...
        zeros(prob_matrix_height, prob_matrix_width, size(f_projections, 2));
    
    % Will contain the eucildean distance betweem the actual projection and
    % the estimated projection.
    dist_matrix = ...
        zeros(prob_matrix_height, prob_matrix_width, size(f_projections, 2));
	
	for k=1:size(f_projections, 2)
		% The current projection.
    	f_proj = f_projections(:, k);

    	% The current orientation for the projection.
	    estimated_orientation = ...
	        Orientation(theta_estimate(k), shift_estimate(k));

	    first_orientation = ...
	    	Orientation(first_estimate_theta(k), first_estimate_shifts(k));

	    % The denominator of the posterior probability.
	    d = denom_posterior_probability_orientation(f_proj,...
            f_image_estimate, estimated_orientation, first_orientation,...
            projection_parameters, noise_estimate, prior_parameters);

	    if d ~= 0
	        for i=-max_angle_err:resolution_angle:max_angle_err
	            for j=-max_shift_err:resolution_space:max_shift_err
	                % The orientation we are currently proposing.
	                given_orientation = ...
	                    Orientation(theta_estimate(k) + i,...
                        shift_estimate(k) + j);

	                % The numerator of the posterior probability.
	                [n, dist] = numer_posterior_probability_orientation(f_proj,...
                        f_image_estimate, first_orientation,...
                        given_orientation, projection_parameters,...
                        noise_estimate, prior_parameters);

	                index_x = round(i/resolution_angle + ...
	                	(max_angle_err + resolution_angle)/resolution_angle);
	                index_y = round(j/resolution_space + ...
	                	(max_shift_err + resolution_space)/resolution_space);

	                % Assign probability to the orientation.
	                prob_matrix(index_x, index_y, k) = n/d;
                    
                    % Assign eucildean distance.
                    dist_matrix(index_x, index_y, k) = dist;
	            end
	        end
	    end
	end
end