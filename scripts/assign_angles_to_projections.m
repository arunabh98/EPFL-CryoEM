function first_theta_estimate = assign_angles_to_projections(...
    f_projections, f_image_estimate, projection_length, output_size)

    first_theta_estimate = zeros(1, size(f_projections, 2));   
    f_image_reshaped = ...
        reshape(f_image_estimate, [output_size, output_size]);
    
    for i=1:size(f_projections, 2)
        min_dist = inf;
        for j=0:1:179
            estimated_projection = project_fourier_alternate(...
				f_image_reshaped, j, 0, projection_length);
            
            if norm(f_projections(:, i) - estimated_projection) < min_dist
                 if sum(ismember(first_theta_estimate(1:i-1), j)) == 0
                    min_dist = norm(f_projections(:, i) - estimated_projection);
                    first_theta_estimate(i) = j;
                 end
            end
        end
    end  
end