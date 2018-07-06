% Get the image.
P = phantom(200);

% Reduce the size of the image for speed.
P = imresize(P, 0.5);

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.05;
max_shift_amplitude = 1;
filename = ...
    '../results/bayesian_estimation/error_angles_and_shifts/server_epfl/5_percent_noise_2/';
num_theta = 180;
max_angle_err = 1;
max_shift_err = 1;
resolution_angle = 0.5;
resolution_space = 1;
no_of_iterations = 3;
mask=ones(size(P));
n = size(P, 1);
L_pad = 260; 

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

% Write the original image.
imwrite(P, strcat(filename, num2str(num_theta), '/original_image.png'));

% Define ground truth angles and take the tomographic projection.
theta = 0:1:179;
[projections, svector] = radon(P, theta);

% The original values.
original_projections = projections;
original_shifts = zeros(size(theta));

% Shift each projection by an unknown amount.
for i=1:size(projections, 2)
    original_shifts(i) = ...
        randi([-max_shift_amplitude, max_shift_amplitude]);
    projections(:, i) = circshift(projections(:, i), original_shifts(i)); 
end
theta_to_write(1, :) = theta;
theta_to_write(7, :) = original_shifts;

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
first_estimate_shifts = original_shifts +...
    randi([-max_shift_err, max_shift_err], 1, num_theta);

% Display the error between the first estimate orientations and the
% correct orientation.
disp(norm(min(abs(first_estimate_theta - theta),...
    180 - abs(first_estimate_theta - theta)), 1));
disp(norm(first_estimate_shifts - original_shifts, 1));

% Begin estimation of the first model.
max_angle_err = 1;
max_shift_err = 1;
prob_matrix_height = (2*max_angle_err)/resolution_angle + 1;
prob_matrix_width = 2*max_shift_err/resolution_space + 1;
prob_matrix = ...
    zeros(prob_matrix_height, prob_matrix_width,...
        size(f_projections, 2)) + 1/(prob_matrix_height*prob_matrix_width);

% Start estimating the image.
fourier_radial = zeros(625, 625);
for i=1:size(prob_matrix, 1)
    for j=1:size(prob_matrix, 2)
        probabilities = squeeze(prob_matrix(i, j, :))';
        prob_f_proj = bsxfun(@mtimes, f_projections, probabilities);

        current_theta = mod(first_estimate_theta  + i*resolution_angle...
            - resolution_angle - max_angle_err, 180);
        current_shift = first_estimate_shifts  + j*resolution_space...
            - resolution_space - max_shift_err;

        fourier_radial = fourier_radial +...
            backproject_fourier_alternate(prob_f_proj, current_theta,...
                current_shift);
    end
end
max_angle_err = 1;
max_shift_err = 1;

f_image_estimate = fourier_radial(:);
first_estimate_model = Ifft2_2_Img(fourier_radial, L_pad);

% Projection and image constants.
output_size = max(size(fourier_radial));
height = size(fourier_radial, 1);
width = size(fourier_radial, 2);
projection_length = size(f_projections, 1);

% Initialize projection parameters object.
projection_parameters = ...
    ProjectionParameters(width, height, output_size, projection_length);

% Show the first estimate image.
figure; imshow(first_estimate_model, []);
disp(norm(first_estimate_model - P));

% Save the first model.
imwrite(first_estimate_model, strcat(filename, num2str(num_theta),...
    '/first_estimate.png'));

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

% Calculate the variance of the prior gaussian distribution.
prior_variance = 10*abs(f_image_estimate).^2;
prior_variance = mean(prior_variance(:));

first_orientation = Orientation(first_estimate_theta, first_estimate_shifts);

% Start the iteration.
theta_estimate = first_estimate_theta;
shift_estimate = first_estimate_shifts;
current_image = first_estimate_model;
error_plot = zeros(no_of_iterations, 1);
weights = ones(size(f_image_estimate));

error_plot(1) = norm(first_estimate_model - P);
for q=1:no_of_iterations
    % The maximum error in angles for this iteration.
    max_angle_err = max(1, max_angle_err - 0.5);

    % Initialize parameters needed for searching in the space.
    prior_parameters = PriorParameters(max_angle_err, max_shift_err,...
        resolution_angle, resolution_space);

    % Calculate probability for each orientation for each projection.
    [prob_matrix, dist_matrix] = calc_prob_for_each_orientation(f_image_estimate,...
        f_projections, theta_estimate, first_estimate_theta,...
        first_estimate_shifts, shift_estimate,...
        noise_estimate, projection_parameters, prior_parameters);

    % Initialize space for correct theta and shifts.
    correct_theta = zeros(size(theta));
    correct_shift = zeros(size(original_shifts));

    % For each projection, choose the right orientation.
    parfor k=1:num_theta
        prob = prob_matrix(:, :, k);
        [~, ind] = max(prob(:));
        [x, y] = ind2sub(size(prob), ind);
        correct_theta(k) = mod(theta_estimate(k)  +...
            x*resolution_angle - resolution_angle - max_angle_err, 180);
        correct_shift(k) = shift_estimate(k) +...
            y*resolution_space - resolution_space - max_shift_err;
    end
    
   % The error between the current orientation and the actual orientations.
    disp(norm(min(abs(correct_theta - theta),...
        180 - abs(correct_theta - theta)), 1));
    disp(norm(correct_shift - original_shifts, 1));

    % Divide the projections with the noise variance.
    denoised_f_projections = ...
        bsxfun(@rdivide, f_projections, (noise_estimate.^2));
    normalization_f_denom = ...
        bsxfun(@rdivide, ones(size(f_projections)), (noise_estimate.^2));

    % Initialize the estimates.
    reconstructed_f_numerator = zeros(size(fourier_radial));
    recons_f_denom = zeros(size(fourier_radial));

    % Start estimating the image.
    for i=1:size(prob_matrix, 1)
        for j=1:size(prob_matrix, 2)
            probabilities = squeeze(prob_matrix(i, j, :))';
            prob_f_proj = bsxfun(@mtimes, denoised_f_projections, probabilities);
            prob_f_noise = bsxfun(@mtimes, normalization_f_denom, probabilities);

            current_theta = mod(theta_estimate  + i*resolution_angle...
                - resolution_angle - max_angle_err, 180);
            current_shift = shift_estimate  + j*resolution_space...
                - resolution_space - max_shift_err;

            reconstructed_f_numerator = reconstructed_f_numerator +...
                backproject_fourier_alternate(prob_f_proj, current_theta,...
                    current_shift);
            recons_f_denom = recons_f_denom +...
                backproject_fourier_alternate(prob_f_noise, current_theta,...
                    current_shift);
        end
    end

    reconstructed_f_numerator = reconstructed_f_numerator(:);
    recons_f_denom = recons_f_denom(:);
    
    % Update the weights on the denominator.
    recons_f_denom_prior = recons_f_denom + 1./prior_variance;
    w = kaiser(size(f_image_estimate, 1), 20);
    recons_f_denom_weighted = ...
        weights./cconv((recons_f_denom_prior.*weights), w, size(f_image_estimate, 1));
    reconstructed_f_image = ...
        reconstructed_f_numerator.*recons_f_denom_weighted*1e4;
    reconstructed_f_image = ...
        reshape(reconstructed_f_image, [output_size, output_size]);

    % Construct the image for this iteration.
    reconstructed_image = Ifft2_2_Img(reconstructed_f_image, L_pad);

    % Display the error for this iteration.
    disp(norm(reconstructed_image - P));

    % Write the image constructed in this reconstruction.
    imwrite(reconstructed_image, strcat(filename, num2str(num_theta),...
        '/reconstructed_image_', num2str(q), '.png'));

    % Note the error and update the parameters for the next iteration.
    error_plot(q + 1) = norm(reconstructed_image - P);
    current_image = reconstructed_image;
    weights = recons_f_denom_weighted;

    % next noise estimate.
    estimated_noise_vector = zeros(projection_length, num_theta);
    parfor k=1:num_theta
        c_proj = project_fourier_alternate(reconstructed_f_image, correct_theta(k),...
            correct_shift(k), projection_length);
        f_proj = f_projections(:, k);
        estimated_noise_vector(:, k) = 0.5*abs(c_proj - f_proj).^2;
    end
    estimated_noise_vector = sqrt(estimated_noise_vector);
    noise_estimate = mean(estimated_noise_vector(:));
    noise_estimate = repmat(noise_estimate, size(f_projections, 1), 1);

    f_image_estimate = reconstructed_f_image(:);
    theta_estimate = correct_theta;
    shift_estimate = correct_shift;

    % Calculate the variance of the prior gaussian distribution.
    prior_variance = 10*abs(f_image_estimate).^2;
    prior_variance = mean(prior_variance(:));
end

% Plot the error as the iteration progresses.
figure; plot(error_plot);
saveas(gcf, strcat(filename, num2str(num_theta), '/error.png'));

% Show the reconstructed image.
figure; imshow(reconstructed_image);

theta_to_write(2, :) = first_estimate_theta;
theta_to_write(3, :) = correct_theta;
theta_to_write(8, :) = first_estimate_shifts;
theta_to_write(9, :) = correct_shift;

% Relative error in theta.
theta_to_write(4, 1) =...
    norm(correct_theta - theta);
theta_to_write(5, 1) =...
    norm(correct_shift - original_shifts);

theta_to_write(6, 1) = norm(current_image - P);

% Write all the parameters and estimated parameters.
csvwrite(strcat(filename, num2str(num_theta),...
    '/thetas.csv'), theta_to_write);

% Write the reconstructed image.
imwrite(reconstructed_image, strcat(filename, num2str(num_theta),...
    '/reconstructed_image.png'));
