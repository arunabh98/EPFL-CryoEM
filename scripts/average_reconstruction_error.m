function noise_estimate = average_reconstruction_error(f_image_estimate,...
    f_projections, theta_estimate, shift_estimate,...
    projection_parameters, prior_parameters, sigmaNoise)

    % Constants
    output_size = projection_parameters.output_size;
    projection_length = projection_parameters.projection_length;
    max_angle_err = prior_parameters.max_angle_err;
    resolution_angle = prior_parameters.resolution_angle;

    fourier_radial = reshape(f_image_estimate, [output_size, output_size]);
    
    U = zeros(size(f_projections, 2), 1);
    U_dist = ...
        zeros(size(f_projections, 2), (2*max_angle_err)/resolution_angle + 1);
    for i=1:size(f_projections, 2)
        f_proj = f_projections(:, i);
        % The orientation specified for this iteration.
        estimated_orientation = Orientation(theta_estimate(i), 0);
        
        % Probability of the projection given the model.
        [U(i), dist_vector] = ...
            calc_prob_projection(f_proj,...
            f_image_estimate, estimated_orientation, sigmaNoise,...
            projection_parameters, prior_parameters);
        U_dist(i, :) = dist_vector;
    end
    
    % Initialize noise estimate.
    noise_estimate = zeros(projection_length, size(f_projections, 2));     
    for j=-max_angle_err:resolution_angle:max_angle_err
        for k=1:size(f_projections, 2)
            % Calculate the projection.
            c_proj = project_fourier_alternate(fourier_radial,...
                theta_estimate(k) + j, shift_estimate(k),...
                projection_length);
            f_proj = f_projections(:, k);
            
            % Orientation of the projection for this iteration.
			orientation = Orientation(theta_estimate(k) + j, 0);
            
            % Calculate the prob of projection given orientation and model.
			[~, dist_proj] = ...
                prob_of_proj_given_orientation(f_proj,...
				f_image_estimate, orientation, sigmaNoise,...
                projection_parameters);
            
            dist_k_orient = U_dist(k, :) - dist_proj;
            prob_k_orient = 1/(sum(exp(vpa(dist_k_orient))));
            
            % Calculate the probabiity of that projection.
            noise_estimate(:, k) = noise_estimate(:, k) +...
                prob_k_orient*0.5*abs(c_proj - f_proj).^2;
        end
    end
    noise_estimate = vpa(mean(noise_estimate(:)));
    noise_estimate = repmat(noise_estimate, projection_length, 1);
end