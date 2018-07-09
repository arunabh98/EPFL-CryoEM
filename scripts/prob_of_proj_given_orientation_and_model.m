function [prob, euler_dist] = prob_of_proj_given_orientation_and_model(f_proj, f_image,...
    orientation, sigmaNoise, projection_parameters)
    projection_length = projection_parameters.projection_length;
    output_size = projection_parameters.output_size;
    f_image = reshape(f_image, [output_size, output_size]);
    c_proj = project_fourier_alternate(f_image,...
        orientation.theta, orientation.shift, projection_length);
    euler_dist = norm(c_proj - f_proj);
    prob = exp(vpa(sum(abs(f_proj - c_proj).^2./(-2*(sigmaNoise.^2))))*1e-3);
end