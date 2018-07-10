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
    euler_dist = norm(c_proj - f_proj);
    cons = vpa(((1/(2*pi))^projection_length)*(1/prod(sigmaNoise/90000))*1e25);
    prob = ...
        cons*exp(vpa(sum((abs(f_proj - c_proj).^2)./(-2*(sigmaNoise))))*1e-2);
end