function [prob, dist_matrix, prior_matrix] = denom_posterior_probability_orientation(...
    f_proj, f_image,...
    estimated_orientation, first_orientation, projection_parameters,... 
    sigmaNoise, prior_parameters)
    
    % Define the constants.
    max_angle_err = prior_parameters.max_angle_err;
    max_shift_err = prior_parameters.max_shift_err;
    resolution_angle = prior_parameters.resolution_angle;
    resolution_space = prior_parameters.resolution_space;
    prob = 0;
    prob_matrix_height = (2*max_angle_err)/resolution_angle + 1;
	prob_matrix_width = 2*max_shift_err/resolution_space + 1;
    
    % Initialize the dist matrix and the prior matrix.
    dist_matrix = ...
        zeros(prob_matrix_height, prob_matrix_width);
    prior_matrix = ...
        zeros(prob_matrix_height, prob_matrix_width);
    
    for i=-max_angle_err:resolution_angle:max_angle_err
        for j=-max_shift_err:resolution_space:max_shift_err
            given_orientation = ...
                Orientation(estimated_orientation.theta + i,...
                estimated_orientation.shift + j);
            % Calculate the probability and the distances.
            [p_a, dist_proj] = prob_of_proj_given_orientation(...
                f_proj, f_image,...
                given_orientation, sigmaNoise, projection_parameters);
            
            % Calculate the probabilities.
            p_b = prior_prob_orientation_given_model(first_orientation,...
                given_orientation, max_angle_err, max_shift_err);
            p = p_b*p_a;
            
            % Calculate the indices in the matrix.
            index_x = round(i/resolution_angle + ...
                (max_angle_err + resolution_angle)/resolution_angle);
            index_y = round(j/resolution_space + ...
                (max_shift_err + resolution_space)/resolution_space);
            
            % Save the distances and the prior probabilities in the matrix.
            prob = prob + p;
            dist_matrix(index_x, index_y) = dist_proj;
            prior_matrix(index_x, index_y) = p_b;
        end
    end
end