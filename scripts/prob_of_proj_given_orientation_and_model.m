function prob = prob_of_proj_given_orientation_and_model(f_proj, f_image,...
    orientation, sigmaNoise, projection_parameters, theta_estimate)
    % S = sliceMatrix(orientation, projection_parameters);
    % c_proj = S*f_image;
    % % a = (1/(2*pi*sigmaNoise^2))^size(f_proj, 1);
    % if size(sigmaNoise, 1) == 1
    % 	sigmaNoise = repmat(sigmaNoise, size(f_proj, 1), 1);
    % end
    projection_length = projection_parameters.projection_length;
    output_size = projection_parameters.output_size;
    f_image = reshape(f_image, [output_size, output_size]);
    c_proj = project_fourier_alternate(f_image,...
        orientation.theta, projection_length);
    prob = exp(vpa(sum(abs(f_proj - c_proj).^2./(-2*(sigmaNoise.^2)))*1e-3));
end