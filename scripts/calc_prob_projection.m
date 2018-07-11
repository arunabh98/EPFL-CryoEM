function [prob, dist_vector] = calc_prob_projection(...
	f_proj, f_image, estimated_orientation, noise_estimate,...
	projection_parameters, prior_parameters)

	% Define prior parameters.
	max_angle_err = prior_parameters.max_angle_err;
	resolution_angle = prior_parameters.resolution_angle;
	angle_estimate = estimated_orientation.theta;
    dist_vector = zeros(1, (2*max_angle_err)/resolution_angle + 1);

	prob = 0;
	for i=-max_angle_err:resolution_angle:max_angle_err
		orientation = ...
	        Orientation(angle_estimate + i, 0);
		[prob_orient, euler_dist] = prob_of_proj_given_orientation(f_proj, f_image,...
            orientation, noise_estimate, projection_parameters);
		prob = prob + vpa(prob_orient);
        
        index_x = round(i/resolution_angle +...
            (max_angle_err + resolution_angle)/resolution_angle);
        dist_vector(index_x) = euler_dist;
	end
end