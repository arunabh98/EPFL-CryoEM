% Get the image.
% P = imread('../images/200px-mickey.jpg');
P = phantom(200);
P = imresize(P, 0.4);
% P = im2double(rgb2gray(P));

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.05;
max_shift_amplitude = 0;
filename = ...
    '../results/bayesian_estimation/unknown_angles_and_shifts/3_percent_noise/';
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);
num_theta = 50;
max_angle_err = 1;
max_shift_err = 0;
resolution_angle = 1;
resolution_space = 1;

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

% Write the original image.
imwrite(P, strcat(filename, num2str(num_theta), '/original_image.png'));

% Define ground truth angles and take the tomographic projection.
theta = datasample(0:179, num_theta);  
[projections, svector] = radon(P,theta);
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
projection_parameters = ProjectionParameters(width, height, output_size, projection_length);

% Initialize parameters needed for searching in the space.
prior_parameters = PriorParameters(max_angle_err, max_shift_err, resolution_angle, resolution_space);

% Add noise to projections.
[projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Transform all entities to the frequency space. 
f_projections = fft(projections);
f_image = fft2(P);
f_image = f_image(:);

% In the first case the angles and shifts will be unknown upto a certain limit.
first_estimate_theta = theta + randi([-max_angle_err, max_angle_err], 1, num_theta);
first_estimate_shifts = original_shifts + randi([-max_shift_err, max_shift_err], 1, num_theta);
first_estimate_model = iradon(projections, first_estimate_theta, output_size);

% Show the first estimate image.
% figure; imshow(first_estimate_model, []);
disp(norm(first_estimate_model - P));

% The initial quantities.
image_estimate = fft2(first_estimate_model);
image_estimate = image_estimate(:);

% The initial noise estimate.
estimated_noise_vector = zeros(projection_length, num_theta);
for k=1:num_theta
    given_orientation = ...
        Orientation(first_estimate_theta(k), first_estimate_shifts(k));
    S = sliceMatrix(given_orientation, projection_parameters);
    c_proj = S*image_estimate;
    f_proj = f_projections(:, k);
    for j=1:projection_length
        estimated_noise_vector(j, k) = 0.5*norm(c_proj(j) - f_proj(j))^2;
    end
end
estimated_noise_vector = sqrt(estimated_noise_vector);
noise_estimate = mean(estimated_noise_vector(:));

% Calculate the variance of the prior gaussian distribution.
prior_variance = zeros(size(image_estimate));
for i=1:size(image_estimate)
    prior_variance(i) = 0.5*norm(image_estimate(i))^2;
end
prior_variance = mean(prior_variance)*4;

% Calculate the noise variance using an alternate formula.
noise_sum = 0;
noise_square_sum = 0;
for k=1:num_theta
    orientation = Orientation(first_estimate_theta(k), first_estimate_shifts(k));
    other_orientation = Orientation(first_estimate_theta(k) - 1, first_estimate_shifts(k));
    other_one = Orientation(first_estimate_theta(k) + 1, first_estimate_shifts(k));
    S_1 = sliceMatrix(orientation, projection_parameters);
    S_2 = sliceMatrix(other_orientation, projection_parameters);
    S_3 = sliceMatrix(other_one, projection_parameters);
    c_proj_1 = S_1*image_estimate;
    c_proj_2 = S_2*image_estimate;
    c_proj_3 = S_3*image_estimate;
    f_proj = f_projections(:, k);
    noise_sum = noise_sum + sum(f_proj - c_proj_1);
    noise_square_sum = noise_square_sum + sum((f_proj - c_proj_1).^2);
    disp(norm(c_proj_1-f_proj));
    disp(norm(c_proj_2-f_proj));
    disp(norm(c_proj_3-f_proj));
end
noise_variance = real((noise_square_sum - (noise_sum/num_theta))/(num_theta - 1));

% for q=1:20

%     % Calculate the probabilities of each orientation for each projection.
%     prob_matrix = zeros((2*max_angle_err)/resolution_angle + 1, 2*max_shift_err/resolution_space + 1, num_theta);

%     for k=1:num_theta
%         % The current projection.
%         f_proj = f_projections(:, k);

%         % The current orientation for the projection.
%         estimated_orientation = ...
%             Orientation(first_estimate_theta(k), first_estimate_shifts(k));

%         % The denominator of the posterior probability.
%         d = denom_posterior_probability_orientation(f_proj, image_estimate,...
%             estimated_orientation, projection_parameters, noise_estimate,...
%             prior_parameters);

%         if d ~= 0
%             for i=-max_angle_err:resolution_angle:max_angle_err
%                 for j=-max_shift_err:resolution_space:max_shift_err
%                     % The orientation we are currently proposing.
%                     given_orientation = ...
%                         Orientation(first_estimate_theta(k) + i, first_estimate_shifts(k) + j);

%                     % The numerator of the posterior probability.
%                     n = numer_posterior_probability_orientation(f_proj, image_estimate,...
%                         estimated_orientation, given_orientation, ...
%                         projection_parameters, noise_estimate, prior_parameters);

%                     prob_matrix(i/resolution_angle + (max_angle_err + resolution_angle)/resolution_angle,...
%                         j/resolution_space + (max_shift_err + resolution_space)/resolution_space, k) = n/d;
%                 end
%             end
%         end
%     end

%     modified_f_projections = f_projections/(noise_estimate^2);
%     numer_image = zeros(size(P));

%     for i=-max_angle_err:resolution_angle:max_angle_err
%         for j=-max_shift_err:resolution_space:max_shift_err
%             prob_vector = prob_matrix(i/resolution_angle + (max_angle_err + resolution_angle)/resolution_angle,...
%                 j/resolution_space + (max_shift_err + resolution_space)/resolution_space, :);

%             prob_vector = squeeze(prob_vector)';
%             prob_vector = repmat(prob_vector, size(modified_f_projections, 1), 1);
            
%             modified_f_projections = modified_f_projections;
%             modified_projections = real(ifft(modified_f_projections)).*prob_vector;

%             for k=1:num_theta
%                 modified_projections(:, k) = circshift(modified_projections(:, k), first_estimate_shifts(k) + j);
%             end

%             numer_image = numer_image +...
%                 iradon(modified_projections, first_estimate_theta(k), output_size)*resolution_angle*resolution_space;
%         end
%     end
% end

% imshow(numer_image, []);

for q=1:20

    % Calculate the probabilities of each orientation for each projection.
    prob_matrix = zeros((2*max_angle_err)/resolution_angle + 1, 2*max_shift_err/resolution_space + 1, num_theta);

    modified_f_projections = f_projections/(noise_estimate^2);
    modified_projections = real(ifft(modified_f_projections));
    numerator = zeros(output_size, output_size);

    normalization_f_noise = ones(size(f_projections))/(noise_estimate^2);
    normalization_noise = real(ifft(normalization_f_noise));
    denom = zeros(output_size, output_size);

    next_noise_vector = zeros(projection_length, num_theta);

    for k=1:num_theta
        % The current projection.
        f_proj = f_projections(:, k);
        
        % The current orientation for the projection.
        estimated_orientation = ...
            Orientation(first_estimate_theta(k), first_estimate_shifts(k));

        % The denominator of the posterior probability.
        d = denom_posterior_probability_orientation(f_proj, image_estimate,...
            estimated_orientation, projection_parameters, noise_estimate,...
            prior_parameters);

        if d ~= 0
            for i=-max_angle_err:resolution_angle:max_angle_err
                for j=-max_shift_err:resolution_space:max_shift_err
                    % The orientation we are currently proposing.
                    given_orientation = ...
                        Orientation(first_estimate_theta(k) + i, first_estimate_shifts(k) + j);

                    % The numerator of the posterior probability.
                    n = numer_posterior_probability_orientation(f_proj, image_estimate,...
                        estimated_orientation, given_orientation, ...
                        projection_parameters, noise_estimate, prior_parameters);

                    prob_matrix(i/resolution_angle + (max_angle_err + resolution_angle)/resolution_angle,...
                        j/resolution_space + (max_shift_err + resolution_space)/resolution_space, k) = n/d;
                    prob = prob_matrix(i/resolution_angle + (max_angle_err + resolution_angle)/resolution_angle,...
                        j/resolution_space + (max_shift_err + resolution_space)/resolution_space, k);

                    current_f_proj = modified_f_projections(:, k);
                    current_f_proj = current_f_proj.*prob;
                    current_proj = real(ifft(current_f_proj));
                    current_proj = circshift(current_proj, first_estimate_shifts(k) + j);
                    numerator = numerator +...
                        iradon([current_proj current_proj], [first_estimate_theta(k) + i first_estimate_theta(k) + i], 'pchip', 'Shepp-Logan', output_size)*resolution_angle*resolution_space;

                    current_proj = normalization_noise(:, k);
                    current_proj = circshift(current_proj, first_estimate_shifts(k) + j);
                    denom = denom +...
                        iradon([current_proj current_proj], [first_estimate_theta(k) + i first_estimate_theta(k) + i], 'pchip', 'Shepp-Logan', output_size)*resolution_angle*resolution_space;

                    S = sliceMatrix(given_orientation, projection_parameters);
                    c_proj = S*image_estimate;
                    for r=1:projection_length
                        next_noise_vector(r, k) = next_noise_vector(r, k) +...
                            0.5*prob*(norm(c_proj(r) - f_proj(r))^2)*resolution_angle*resolution_space;
                    end
                end
            end
        end
    end 

    numerator = fft2(numerator);
    numerator = numerator(:);

    denom = fft2(denom);
    denom = denom(:);

    next_image_estimate = numerator./(denom + prior_variance);
    img = reshape(next_image_estimate, [output_size, output_size]);
    img = real(ifft2(img));
    disp(norm(img-P));
    figure; imshow(img, []);

    % Estimates for the nest iteration.
    next_noise_vector = sqrt(next_noise_vector);
    noise_estimate = mean(next_noise_vector(:));
    image_estimate = next_image_estimate;

    % Calculate the variance of the prior gaussian distribution.
    prior_variance = zeros(size(image_estimate));
    for i=1:size(image_estimate)
        prior_variance(i) = 0.5*norm(image_estimate(i))^2;
    end
    prior_variance = mean(prior_variance)*4;
end