% Get the image.
P = phantom(200);

% Pad the image with a fixed boundary of 3 pixels.
% P = padarray(P, [3th, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.02;
max_shift_amplitude = 0;
filename = ...
    '../results/bayesian_estimation/error_angles_and_shifts/5_percent_noise/';
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);
num_theta = 180;
max_angle_err = 1;
max_shift_err = 0;
resolution_angle = 1;
resolution_space = 1;
no_of_iterations = 1;

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

% Write the original image.
imwrite(P, strcat(filename, num2str(num_theta), '/original_image.png'));

% Define ground truth angles and take the tomographic projection.
% theta = datasample(0:0.5:359.5, num_theta);  
theta = 0:179;
[projections, svector] = radon(P, theta);
original_projections = projections;
original_shifts = zeros(size(theta));

% Shift each projection by an unknown amount.
for i=1:size(projections, 2)
    original_shifts(i) = ...
        randi([-max_shift_amplitude, max_shift_amplitude]);
    projections(:, i) = circshift(projections(:, i), original_shifts(i)); 
end
theta_to_write(1, :) = theta;
theta_to_write(6, :) = original_shifts;

% The length of the projections.
projection_length = size(projections, 1);

% Initialize projection parameters object.
projection_parameters = ...
    ProjectionParameters(width, height, output_size, projection_length);

% Initialize parameters needed for searching in the space.
prior_parameters = PriorParameters(max_angle_err, max_shift_err,...
    resolution_angle, resolution_space);

% Add noise to projections.
[projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Transform all entities to the frequency space. 
f_projections = fftshift(fft(ifftshift(projections)));
f_image = fftshift(fft2(ifftshift(P)));
f_image_vector = f_image(:);

% In the first case the angles and shifts will be unknown upto a 
% certain limit.
first_estimate_theta = theta +...
    randi([-max_angle_err, max_angle_err], 1, num_theta);
first_estimate_shifts = original_shifts +...
    randi([-max_shift_err, max_shift_err], 1, num_theta);
first_estimate_model = ...
    iradon(projections, first_estimate_theta, output_size);

% [fourier_radial, first_estimate_model_alternate] = ...
%     Direct_Fourier_Reconst(P, first_estimate_theta, projections);
first_estimate_model_alternate = flipdim(first_estimate_model_alternate, 1);
slope = max(first_estimate_model(:)) - min(first_estimate_model(:));
c = min(first_estimate_model(:));

% Show the first estimate image.
% figure; imshow(first_estimate_model, []);
% figure; imshow(first_estimate_model_alternate, []);
disp(norm(first_estimate_model - P));
% disp(norm(first_estimate_model_alternate - P));

% The initial quantities.
f_image_estimate = fftshift(fft2(ifftshift(first_estimate_model)));
f_image_estimate = f_image_estimate(:);

% f_image_estimate_alter = fftshift(fft2(ifftshift(first_estimate_model_alternate)));
% f_image_estimate_alter = f_image_estimate_alter(:);

% Just a necessary check.
original_recons_image = reshape(f_image_vector, [output_size, output_size]);
original_recons_image = ifftshift(ifft2(fftshift(original_recons_image)));
figure; imshow(original_recons_image, []);

direct_recons_image = reshape(f_image_estimate, [output_size, output_size]);
direct_recons_image = ifftshift(ifft2(fftshift(direct_recons_image)));
figure; imshow(direct_recons_image, []);

% slice_recons_image = reshape(f_image_estimate_alter, [output_size, output_size]);
% slice_recons_image = ifftshift(ifft2(fftshift(slice_recons_image)));
% figure; imshow(slice_recons_image, []);

% Initial noise_estimate.
noise_estimate = repmat(1/sigmaNoise, [projection_length, 1]);

% Actual noise estimate.
actual_noise_estimate = zeros(projection_length, num_theta);
for k=1:num_theta
    given_orientation = ...
        Orientation(theta(k), original_shifts(k));
    S = sliceMatrix(given_orientation, projection_parameters);
    c_proj = S*f_image_vector;
    f_proj = f_projections(:, k);
    actual_noise_estimate(:, k) = 0.5*abs(c_proj - f_proj).^2;
end
actual_noise_estimate = sqrt(actual_noise_estimate);
actual_noise_estimate = mean(actual_noise_estimate, 2);

% Actual prior estimate.
actual_prior_variance = 2*abs(f_image_vector).^2;

% Calculate the variance of the prior gaussian distribution.
prior_variance = 2*abs(f_image_estimate).^2;
% prior_variance_alter = 2*abs(f_image_estimate_alter).^2;

% Start the iteration.
theta_estimate = first_estimate_theta;
shift_estimate = first_estimate_shifts;
current_image = first_estimate_model;

error_plot = zeros(no_of_iterations, 1);
weights = ones(size(f_image_estimate));

% The actual values.
noise_estimate = actual_noise_estimate;
f_image_estimate = f_image_vector;
prior_variance = actual_prior_variance;

for q=1:no_of_iterations
    % Calculate probability for each orientation for each projection.
    prob_matrix = calc_prob_for_each_orientation(f_image_estimate,...
        f_projections, theta_estimate, shift_estimate, noise_estimate,...
        projection_parameters, prior_parameters);

    % Initialize space for correct theta and shifts.
    correct_theta = zeros(size(theta));
    correct_shift = zeros(size(original_shifts));

    % For each projection, choose the right orientation.
    for k=1:num_theta
        prob = prob_matrix(:, :, k);
        [~, ind] = max(prob(:));
        [x, y] = ind2sub(size(prob), ind);
        correct_theta(k) = theta_estimate(k)  +...
            x*resolution_angle - resolution_angle - max_angle_err;
        correct_shift(k) = shift_estimate(k) +...
            y*resolution_space - resolution_space - max_shift_err;
    end

    % Start estimating the image.
    denoised_f_projections = ...
        bsxfun(@rdivide, f_projections, (noise_estimate.^2));
    normalization_f_denom = ...
        bsxfun(@rdivide, ones(size(f_projections)), (noise_estimate.^2));

    denoised_projections = real(ifftshift(ifft(fftshift(denoised_f_projections))));
    normalization_denom = real(ifftshift(ifft(fftshift(normalization_f_denom))));

    % Correct the shift in the projections.
    corrected_projections = zeros(size(projections));
    normalization_denom_image = zeros(size(projections));
    for i=1:num_theta
        corrected_projections(:, i) = ...
            circshift(denoised_projections(:, i), -correct_shift(i));
        normalization_denom_image(:, i) = ...
            circshift(normalization_denom(:, i), -correct_shift(i));
    end

    % Reconstruct the numerator and denominator images.
    reconstructed_numerator =...
        iradon(corrected_projections, correct_theta, output_size);
    reconstructed_denominator =...
        iradon(normalization_denom_image, correct_theta, output_size);

    reconstructed_f_numerator = fftshift(fft2(ifftshift(reconstructed_numerator)));
    reconstructed_f_numerator = reconstructed_f_numerator(:);
    recons_f_denom = fftshift(fft2(ifftshift(reconstructed_denominator)));
    recons_f_denom = recons_f_denom(:);

    error = inf;
    % for j=1e-6:1e-5:1e-4
        for i=-50:1:50
            recons_f_denom_prior = recons_f_denom + (1./prior_variance);
            w = kaiser(size(f_image_estimate, 1), i);
            recons_f_denom_weighted = ...
                weights./conv((recons_f_denom_prior.*weights), w, 'same');
            reconstructed_f_image = ...
                reconstructed_f_numerator.*recons_f_denom_weighted;
            reconstructed_f_image = ...
                reshape(reconstructed_f_image, [output_size, output_size]);

            estimated_image = mat2gray(real(ifftshift(ifft2(fftshift(reconstructed_f_image)))));
            estimated_image = slope.*estimated_image + c;

            if norm(estimated_image - current_image) < error
                error = norm(estimated_image - current_image);
                reconstructed_image = estimated_image;
                % optimal_variance = j;
                optimal_b = i;
                perfect_weights = recons_f_denom_weighted;
            end
        end
    % end
    disp(norm(reconstructed_image - P));
    error_plot(q) = norm(reconstructed_image - P)/norm(P);
    % disp(optimal_variance);
    disp(optimal_b);
    current_image = reconstructed_image;
    reconstructed_f_image = fftshift(fft2(ifftshift(reconstructed_image)));
    weights = recons_f_denom_weighted;
    figure; imshow(reconstructed_image, []);

    % next noise estimate.
    estimated_noise_vector = zeros(projection_length, num_theta);
    for k=1:num_theta
        given_orientation = ...
            Orientation(theta_estimate(k), shift_estimate(k));
        S = sliceMatrix(given_orientation, projection_parameters);
        c_proj = S*f_image_estimate;
        f_proj = f_projections(:, k);
        estimated_noise_vector(:, k) = 0.5*abs(c_proj - f_proj).^2;
    end
    estimated_noise_vector = sqrt(estimated_noise_vector);

    f_image_estimate = reconstructed_f_image(:);
    theta_estimate = correct_theta;
    shift_estimate = correct_shift;
    % noise_estimate = mean(estimated_noise_vector, 2);

    % Calculate the variance of the prior gaussian distribution.
    % prior_variance = 2*abs(f_image_estimate).^2;
end

% Plot the error as the iteration progresses.
figure; plot(error_plot);
saveas(gcf, strcat(filename, num2str(num_theta), '/error.png'));

% Show the reconstructed image.
figure; imshow(reconstructed_image, []);

theta_to_write(2, :) = correct_theta;
theta_to_write(7, :) = correct_shift;

% Relative error in theta.
theta_to_write(3, 1) =...
    norm(correct_theta - theta)/norm(theta);
theta_to_write(4, 1) =...
    norm(correct_shift - original_shifts)/norm(original_shifts);

theta_to_write(5, 1) = norm(current_image - P)/norm(P);

% Write all the parameters and estimated parameters.
csvwrite(strcat(filename, num2str(num_theta),...
    '/thetas.csv'), theta_to_write);

% Write the reconstructed image.
imwrite(reconstructed_image, strcat(filename, num2str(num_theta),...
    '/reconstructed_image.png'));


% error = inf;
% optimal_b = 0;
% % optimal_variance = 0;
% % for j=1e-5:1e-5:1e-3
%     for i=-50:0.5:50
%         recons_f_denom = recons_f_denom + 2e-5;
%         w = kaiser(size(f_image_estimate, 1), i);

%         recons_f_denom = S./cconv((recons_f_denom.*S), w, size(f_image_estimate, 1));

%         reconstructed_f_image = reconstructed_f_numerator./recons_f_denom;
%         reconstructed_f_image = reshape(reconstructed_f_image, [output_size, output_size]);
%         reconstructed_image = real(ifft2(reconstructed_f_image));

%         if norm(reconstructed_image - P) < error
%             error = norm(reconstructed_image - P);
%             calc_image = reconstructed_image;
%             optimal_b = i;
%             optimal_variance = j;
%         end
%     end
% % end;
% disp(error);
% disp(optimal_b);
% % disp(optimal_variance);
% figure; imshow(calc_image);
% % % next noise estimate.
% % estimated_noise_vector = zeros(projection_length, num_theta);
% % for k=1:num_theta
% %     given_orientation = ...
% %         Orientation(theta(k), original_shifts(k));
% %     S = sliceMatrix(given_orientation, projection_parameters);
% %     c_proj = S*f_image_vector;
% %     f_proj = f_projections(:, k);
% %     estimated_noise_vector(:, k) = 0.5*abs(c_proj - f_proj).^2;
% % end
% % estimated_noise_vector = sqrt(estimated_noise_vector);
% % close all;
