digits(100);
% Get the image.
P = phantom(200);

% Reduce the size of the image for speed.
P = imresize(P, 0.5);

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.05;
filename = ...
	'../results/cryoSPARC/5_percent_noise/';
num_theta = 180;
max_shift_err = 0;
max_angle_err = 30;
resolution_angle = 1;
resolution_space = 1;
L_pad = 260; 
no_of_iterations = 10;
momentum_parameter = 0.9;
gamma = 0.9999;
number_of_samples = zeros(no_of_iterations, 1);
number_of_samples(1) = 30;

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

% Write the original image.
imwrite(P, strcat(filename, num2str(num_theta), '/original_image.png'));

% Define ground truth angles and take the tomographic projection.
theta = 0:1:179;
[projections, svector] = radon(P, theta);

% The original values.
original_projections = projections;
theta_to_write(1, :) = theta;

% Add noise to projections.
[projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Transform all entities to the frequency space.
f_projections = ifftshift(projections,1);
f_projections = fft(f_projections, [ ], 1);
f_projections = fftshift(f_projections, 1); % put DC central after filtering 

% In the first case the angles and shifts will be unknown upto a 
% certain limit.
first_estimate_theta = mod(theta +...
	randi([-max_angle_err, max_angle_err], 1, num_theta), 180);
first_estimate_shifts = zeros(size(theta));

% Size of fourier domain.
f_size = 625;

% Initialize parameters needed for searching in the space.
prior_parameters = PriorParameters(max_angle_err, max_shift_err,...
    resolution_angle, resolution_space);

% Randomly initialize the initial model with std_deviation = 5
std_deviation = 5;
f_image_estimate = ...
	std_deviation*randn(f_size*f_size, 1) + i*std_deviation*randn(f_size*f_size, 1);
prior_variance = repmat(std_deviation.^2*50, size(f_image_estimate, 1), 1);

% Make the first estimate in the fourier and the spatial domain.
fourier_radial = reshape(f_image_estimate, [f_size f_size]);
first_estimate_model = Ifft2_2_Img(fourier_radial, L_pad);

% Write the first estimate image.
imwrite(first_estimate_model,...
    strcat(filename, num2str(num_theta), '/original_image.png'));

% Projection and image constants.
output_size = max(size(fourier_radial));
height = size(fourier_radial, 1);
width = size(fourier_radial, 2);
projection_length = size(f_projections, 1);

first_noise_estimate = repmat(1/(sigmaNoise.^2)*3100, projection_length, 1);

% Initialize projection parameters object.
projection_parameters = ...
	ProjectionParameters(width, height, output_size, projection_length);

% Noise estimate.
noise_estimate = zeros(projection_length, 1);
noise_estimate(:, 1) = average_reconstruction_error(f_image_estimate,...
    f_projections, first_estimate_theta, first_estimate_shifts,...
    projection_parameters, prior_parameters, first_noise_estimate);

% Show the first estimate image.
figure; imshow(first_estimate_model, []);
disp(norm(first_estimate_model - P));

% Save the first error.
error_plot = zeros(no_of_iterations + 1, 1);
error_plot(1) = norm(first_estimate_model - P);

% Save the first model.
imwrite(first_estimate_model, strcat(filename, num2str(num_theta),...
	'/first_estimate_image.png'));

% The gradient term for each iteration.
momentum_gradient = zeros(size(f_image_estimate));

% Start the Stochastic Gradient Descent.
for q=1:no_of_iterations
	% Noise estimate for this iteration.
	w_q = zeros(size(noise_estimate, 1), 1);
	sigma_bar_q = zeros(size(noise_estimate, 1), 1);
	w_hat_q = 1*(gamma^q);
	sigma_bar = 1/(sigmaNoise.^2)*3100;
	sigmaHat = 1*sigma_bar;
	parfor i=1:q
		w_q = w_q +...
            ones(size(noise_estimate, 1), 1).*gamma^(q-i)*number_of_samples(i);
		sigma_bar_q = sigma_bar_q + noise_estimate(:, i).*gamma^(q-i);
	end
	noise_variance = ...
		vpa((w_q.*sigma_bar_q + 50*sigma_bar + sigmaHat*w_hat_q)./...
        (w_q + 50 + w_hat_q));
	
	% Randomly select some projections.
	[sel_projections, idx] =...
        datasample(f_projections, number_of_samples(q), 2,...
        'Replace', false);
	selected_angles = first_estimate_theta(idx);
	shifts = zeros(size(selected_angles));

	% The gradient vector for this iteration.
	gradient_vector_iter = zeros(size(f_image_estimate));	
	noisy_step_vector_iter = zeros(size(f_image_estimate));
    
    % Calculate the probability of each projection.
    U = zeros(size(sel_projections, 2), 1);
    U_dist = ...
        zeros(size(sel_projections, 2),...
        (2*max_angle_err)/resolution_angle + 1);
    parfor i=1:size(sel_projections, 2)
        % The orientation specified for this iteration.
        estimated_orientation = Orientation(selected_angles(i), 0);
        
        % Probability of the projection given the model.
        [U(i), dist_vector] = ...
            calc_prob_projection(sel_projections(:, i),...
            f_image_estimate, estimated_orientation, noise_variance,...
            projection_parameters, prior_parameters);
        U_dist(i, :) = dist_vector;
    end

	parfor j=-max_angle_err:resolution_angle:max_angle_err
        
		thetas_iter = selected_angles  + j;
		projections_iter = zeros(size(sel_projections));
		noisy_iter = zeros(size(sel_projections));
        
		for i=1:size(sel_projections, 2)
			% Orientation of the projection for this iteration.
			orientation = Orientation(thetas_iter(i), 0);

			% Calculate the prob of projection given orientation and model.
			[prob_proj, dist_proj] = ...
                prob_of_proj_given_orientation(sel_projections(:, i),...
				f_image_estimate, orientation, noise_variance,...
                projection_parameters);

			% Difference between given projection and calculated projection.
			f_image_reshaped =...
				reshape(f_image_estimate, [output_size, output_size]);
			diff_proj =  ...
                sel_projections(:, i) - project_fourier_alternate(...
				f_image_reshaped, thetas_iter(i), 0, projection_length);
            
            % The probability of each orientation.
			prob_phi = 1/((2*max_angle_err)/resolution_angle + 1);
            
            dist_i_orient = U_dist(i, :) - dist_proj;
            prob_i_orient = 1/(sum(exp(vpa(dist_i_orient))));
            
            % Save projections to find gradient and step-size.
			projections_iter(:, i) = ...
				prob_i_orient.*(diff_proj./noise_variance);	
			noisy_iter(:, i) = ...
				prob_i_orient.*(ones(size(diff_proj))./noise_variance);
		end
		
		% Calculate the gradient due to current estimate.
		gradient_iter = ...
			backproject_fourier_alternate(...
            projections_iter, thetas_iter, shifts);
		gradient_iter = gradient_iter(:);
		
		% Calculate the step size form noise.
		noisy_step_size = ...
			backproject_fourier_alternate(noisy_iter, thetas_iter, shifts);
		noisy_step_size = noisy_step_size(:);

		% Compute the total gradient.
		gradient_vector_iter = gradient_vector_iter + gradient_iter;	
		noisy_step_vector_iter = noisy_step_vector_iter + noisy_step_size;	
	end
	
	% Calculate the gradient due to the prior distribution.
	prior_gradient = f_image_estimate./(-prior_variance);
	
	% Calculate the total gradient for this iteration.
	gradient_vector_iter = ...
		gradient_vector_iter*size(projections, 2)/size(sel_projections, 2) +...
		prior_gradient;
	
	% Calculate the step size.
	step_size = 1/max(noisy_step_vector_iter);
	
	% Gradient with momentum.
	momentum_gradient = ...
		momentum_parameter*momentum_gradient +...
		(1 - momentum_parameter)*step_size*gradient_vector_iter;
	
	% Calculate the next estimate.
	f_image_estimate = f_image_estimate + momentum_gradient;
	f_image_estimate_reshaped = ...
		reshape(f_image_estimate, [output_size, output_size]);
	
	% next noise estimate.
	next_noise_estimate = average_reconstruction_error(f_image_estimate,...
    f_projections, first_estimate_theta, first_estimate_shifts,...
    projection_parameters, prior_parameters, noise_variance);
    noise_estimate = [noise_estimate next_noise_estimate];
	
	% Display the image and calculate the error.
	image_estimate = Ifft2_2_Img(f_image_estimate_reshaped, L_pad);
	figure; imshow(image_estimate, []);
	disp(norm(image_estimate - P));
    
    % Update the number of samples for the next iteration.
    if norm(image_estimate - P) < 7
        number_of_samples(q+1) = 100;
    else
        number_of_samples(q+1) = 30;
    end
    
    % Save the model and the image in the current iteration.
    error_plot(q+1) = norm(image_estimate - P);
    imwrite(image_estimate, strcat(filename, num2str(num_theta),...
        '/reconstructed_image_', num2str(q), '.png'));
end

% Save the error plot.
figure; plot(error_plot);
saveas(gcf, strcat(filename, num2str(num_theta), '/error.png'));

