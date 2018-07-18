digits(100);
% Get the image.
P = phantom(200);

% Reduce the size of the image for speed.
P = imresize(P, 0.5);

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants
filename = ...
    '../results/refine_image/5_percent_noise/';
num_theta = 180;
max_angle_err = 20;
max_shift_err = 0;
resolution_angle = 4;
resolution_space = 1;
no_of_iterations = 10;
L_pad = 260; 
output_size = 625;
num_projections = num_theta;
theta = 0:1:179;

noisy_theta = mod(theta +...
    randi([-max_angle_err, max_angle_err], 1, num_theta), 180);

% Load the projections and initial estimate of image.
saved_filename = ...
    '../results/cryoSPARC/5_percent_noise/office/';
saved_projections = ...
    matfile(strcat(saved_filename, num2str(num_projections),...
    '/f_projections.mat'));
saved_image_estimate = ...
    matfile(strcat(saved_filename, num2str(num_projections),...
    '/f_image_estimate.mat'));
saved_noise_estimate = ...
    matfile(strcat(saved_filename, num2str(num_projections),...
    '/noise_estimate.mat'));

f_projections = saved_projections.f_projections;
f_image_estimate = saved_image_estimate.f_image_estimate;
noise_estimate = ones(size(f_projections, 1), 1);
fourier_radial = reshape(f_image_estimate, [output_size, output_size]);

% Low pass the image in the fourier domain.
[N, M] = size(fourier_radial);
%Sampling intervals 
dx = 1; 
dy = 1; 
%Characteristic wavelengths 
KX0 = (mod(1/2 + (0:(M-1))/M, 1) - 1/2); 
KX1 = KX0 * (2*pi/dx); 
KY0 = (mod(1/2 + (0:(N-1))/N, 1) - 1/2); 
KY1 = KY0 * (2*pi/dx); 
[KX,KY] = meshgrid(KX1,KY1); 
%Filter formulation 
K0 = 3.5;
lpf = ~(KX.*KX + KY.*KY < K0^2);
fourier_radial = fourier_radial.*lpf;
f_image_estimate = fourier_radial(:);

% The first image.
first_estimate_model = Ifft2_2_Img(fourier_radial, L_pad);

% Show the first estimate image.
figure; imshow(first_estimate_model, []);
disp(norm(first_estimate_model - P));

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

% In the first case the angles and shifts will be unknown upto a 
% certain limit.
first_estimate_theta = assign_angles_to_projections(...
    f_projections, f_image_estimate, size(f_projections, 1), output_size,...
    max_angle_err, resolution_angle, noisy_theta);

% Display the error between the first estimate orientations and the
% correct orientation.
disp(norm(min(abs(first_estimate_theta - theta),...
    180 - abs(first_estimate_theta - theta)), 1));

% Zero shift
first_estimate_shifts = zeros(1, num_theta);
original_shifts = zeros(1, num_theta);

% Begin estimation of the first model.
prob_matrix_height = (2*max_angle_err)/resolution_angle + 1;
prob_matrix_width = 2*max_shift_err/resolution_space + 1;
prob_matrix = ...
    zeros(prob_matrix_height, prob_matrix_width,...
        size(f_projections, 2)) + 1/(prob_matrix_height*prob_matrix_width);

% Start estimating the image.
fourier_radial = zeros(625, 625);
parfor i=1:size(prob_matrix, 1)
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

f_image_estimate = fourier_radial(:);
first_estimate_model = Ifft2_2_Img(fourier_radial, L_pad);

% Show the first estimate image.
figure; imshow(first_estimate_model, []);
disp(norm(first_estimate_model - P));

% Projection and image constants.
output_size = max(size(fourier_radial));
height = size(fourier_radial, 1);
width = size(fourier_radial, 2);
projection_length = size(f_projections, 1);

% Initialize projection parameters object.
projection_parameters = ...
    ProjectionParameters(width, height, output_size, projection_length);

% Save the first model.
imwrite(first_estimate_model, strcat(filename, num2str(num_theta),...
    '/first_estimate.png'));

% Calculate the variance of the prior gaussian distribution.
prior_variance = 20*abs(f_image_estimate).^2;
prior_variance = mean(prior_variance(:));

% Start the iteration.
theta_estimate = first_estimate_theta;
shift_estimate = first_estimate_shifts;
current_image = first_estimate_model;
error_plot = zeros(no_of_iterations, 1);
weights = ones(size(f_image_estimate));

error_plot(1) = norm(first_estimate_model - P);
for q=1:no_of_iterations
    % The maximum error in angles for this iteration.
    max_angle_err = max(1, max_angle_err - 5);
    resolution_angle = max(1, resolution_angle - 1);

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
        bsxfun(@rdivide, f_projections, double(noise_estimate));
    normalization_f_denom = ...
        bsxfun(@rdivide, ones(size(f_projections)), double(noise_estimate));

    % Initialize the estimates.
    reconstructed_f_numerator = zeros(size(fourier_radial));
    recons_f_denom = zeros(size(fourier_radial));

    % Start estimating the image.
    parfor i=1:size(prob_matrix, 1)
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
    reconstructed_f_image = reconstructed_f_numerator;
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
    noise_estimate = ones(size(f_projections, 1), 1);
    
    % Next estimate.
    f_image_estimate = reconstructed_f_image(:);
    theta_estimate = correct_theta;
    shift_estimate = correct_shift;

    % Calculate the variance of the prior gaussian distribution.
    prior_variance = 20*abs(f_image_estimate).^2;
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
