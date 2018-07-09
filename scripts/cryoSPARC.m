% Get the image.
P = phantom(200);

% Reduce the size of the image for speed.
P = imresize(P, 0.5);

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.05;
filename = ...
	'../results/bayesian_estimation/error_angles_and_shifts/server_epfl/5_percent_noise/';
num_theta = 180;
max_shift_err = 0;
max_angle_err = 1;
resolution_angle = 1;
resolution_space = 1;
L_pad = 260; 
no_of_iterations = 10;
momentum_parameter = 0.9;

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

% Randomly initialize the initial model with std_deviation = 5
std_deviation = 5;
f_image_estimate = ...
	std_deviation*randn(f_size*f_size, 1) + i*std_deviation*randn(f_size*f_size, 1);
prior_variance = repmat(std_deviation.^2, size(f_image_estimate, 1), 1);

% Make the first estimate in the fourier and the spatial domain.
fourier_radial = reshape(f_image_estimate, [f_size f_size]);
first_estimate_model = Ifft2_2_Img(fourier_radial, L_pad);

% Projection and image constants.
output_size = max(size(fourier_radial));
height = size(fourier_radial, 1);
width = size(fourier_radial, 2);
projection_length = size(f_projections, 1);

% Noise estimate.
noise_estimate = zeros(projection_length, num_theta);
parfor k=1:num_theta
	c_proj = project_fourier_alternate(fourier_radial,...
		first_estimate_theta(k), first_estimate_shifts(k), projection_length);
	f_proj = f_projections(:, k);
	noise_estimate(:, k) = 0.5*abs(c_proj - f_proj).^2;
end
noise_estimate = sqrt(noise_estimate);
noise_estimate = mean(noise_estimate(:));
noise_estimate = repmat(noise_estimate, size(f_projections, 1), 1);

% Initialize projection parameters object.
projection_parameters = ...
	ProjectionParameters(width, height, output_size, projection_length);

% Show the first estimate image.
figure; imshow(first_estimate_model, []);
disp(norm(first_estimate_model - P));

% Save the first model.
imwrite(first_estimate_model, strcat(filename, num2str(num_theta),...
	'/first_estimate.png'));

% The gradient term for each iteration.
momentum_gradient = zeros(size(f_image_estimate));

% Start the Stochastic Gradient Descent.
for q=1:no_of_iterations
	% Randomly select some projections.
	[sel_projections, idx] = datasample(f_projections, 30, 2, 'Replace', false);
	selected_angles = first_estimate_theta(idx);
	shifts = zeros(size(selected_angles));

	% Initialize parameters needed for searching in the space.
	prior_parameters = PriorParameters(max_angle_err, max_shift_err,...
		resolution_angle, resolution_space);

	% The gradient vector for this iteration.
	gradient_vector_iter = zeros(size(f_image_estimate));
    
    noisy_step_vector_iter = zeros(size(f_image_estimate));

	for j=-max_angle_err:resolution_angle:max_angle_err
		thetas_iter = mod(selected_angles  + j, 180);

		projections_iter = zeros(size(sel_projections));
        noisy_iter = zeros(size(sel_projections));
		for i=1:size(sel_projections, 2)
			% Orientation of the projection for this iteration.
			orientation = Orientation(thetas_iter(i), 0);

			% The orientation specified for this iteration.
			estimated_orientation = Orientation(selected_angles(i), 0);

			% Probability of the projection given the model.
			U_i = calc_prob_projection(sel_projections(:, i), f_image_estimate,...
				estimated_orientation, noise_estimate, projection_parameters,...
				prior_parameters);

			% Calculate the prob of projection given orientation and model.
			prob_proj = prob_of_proj_given_orientation_and_model(sel_projections(:, i),...
				f_image_estimate, orientation, noise_estimate, projection_parameters);

			% Difference between given projection and calculated projection.
			f_image_reshaped =...
				reshape(f_image_estimate, [output_size, output_size]);
			diff_proj =  sel_projections(:, i) - project_fourier_alternate(...
				f_image_reshaped, thetas_iter(i), 0, projection_length);

			prob_phi = (1/((2*max_angle_err)/resolution_angle + 1));
            projections_iter(:, i) = ...
				(prob_proj/U_i).*(diff_proj./noise_estimate).*prob_phi;
		
            noisy_iter(:, i) = ...
                (prob_proj/U_i).*(ones(size(diff_proj))./noise_estimate).*prob_phi;
        end
		
		% Calculate the gradient due to current estimate.
		gradient_iter = ...
            backproject_fourier_alternate(projections_iter, thetas_iter, shifts);
		gradient_iter = gradient_iter(:);
        
        % Calculate the step size form noise.
        noisy_step_size = ...
            backproject_fourier_alternate(noisy_iter, thetas_iter, shifts);

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
    
    % Display the image and calculate the error.
    image_estimate = Ifft2_2_Img(f_image_estimate_reshaped, L_pad);
    figure; imshow(image_estimate, []);
    disp(norm(image_estimate - P));
end


