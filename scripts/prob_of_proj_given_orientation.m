function [prob, euler_dist] = prob_of_proj_given_orientation(f_proj, f_image,...
    orientation, sigmaNoise, projection_parameters)
    
    % Constants
    projection_length = projection_parameters.projection_length;
    output_size = projection_parameters.output_size;
    
    % Calculate the projection.
    f_image = reshape(f_image, [output_size, output_size]);
    c_proj = project_fourier_alternate(f_image,...
        orientation.theta, orientation.shift, projection_length);
    
    % Calculate the distance and the probability.
    euler_dist = ...
        vpa(sum((abs(f_proj - c_proj).^2)./(-2*(sigmaNoise))));
    prob = ...
        exp(vpa(sum((abs(f_proj - c_proj).^2)./(-2*(sigmaNoise))))*1e-2);
end