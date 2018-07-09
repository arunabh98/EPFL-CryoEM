function prob = calc_prob_projection(...
	f_proj, f_image, estimated_orientation, noise_estimate,...
	projection_parameters, prior_parameters)

	% Define prior parameters.
	max_angle_err = prior_parameters.max_angle_err;
	resolution_angle = prior_parameters.resolution_angle;
	angle_estimate = estimated_orientation.theta;

	prob = 0;
	for i=-max_angle_err:resolution_angle:max_angle_err
		orientation = ...
	        Orientation(angle_estimate + i, 0);
		prob_orient = ...
			(1/((2*max_angle_err)/resolution_angle + 1))*...
            prob_of_proj_given_orientation_and_model(f_proj, f_image,...
            orientation, noise_estimate, projection_parameters);
		prob = prob + prob_orient.*(1/((2*max_angle_err)/resolution_angle + 1));
	end
end