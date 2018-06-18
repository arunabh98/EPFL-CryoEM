function prob = denom_posterior_probability_orientation(f_proj, f_image,...
    estimated_orientation, projection_parameters, sigmaNoise, prior_parameters,...
    theta_estimate)
    
    
    max_angle_err = prior_parameters.max_angle_err;
    max_shift_err = prior_parameters.max_shift_err;
    resolution_angle = prior_parameters.resolution_angle;
    resolution_space = prior_parameters.resolution_space;
    prob = 0;
    
    for i=-max_angle_err:resolution_angle:max_angle_err
        for j=-max_shift_err:resolution_space:max_shift_err
            given_orientation = ...
                Orientation(estimated_orientation.theta + i,...
                estimated_orientation.shift + j); 
            p_a = prob_of_proj_given_orientation_and_model(f_proj, f_image,...
                given_orientation, sigmaNoise, projection_parameters, theta_estimate);
            p_b = prior_prob_orientation_given_model(estimated_orientation,...
                given_orientation, max_angle_err, max_shift_err);
            p = p_b*p_a;
            prob = prob + p*resolution_angle*resolution_space;
        end
    end
end