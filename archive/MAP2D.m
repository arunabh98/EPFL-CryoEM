% Get the image.
% P = imread('../images/200px-mickey.jpg');
P = phantom(200);
P = imresize(P, 0.4);
% P = im2double(rgb2gray(P));

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.00;
max_shift_amplitude = 0;
filename = ...
    '../results/bayesian_estimation/unknown_angles_and_shifts/3_percent_noise/';
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);
num_theta = 80;
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
figure; imshow(first_estimate_model, [])

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
prior_variance = prior_variance*4;

for q=1:3
    % The next noise estimate.
    next_noise_vector = zeros(projection_length, num_theta);

    % The next image estimate.
    numerator_image = zeros(output_size*output_size, 1);
    denominator_image = zeros(output_size*output_size, 1);
    % Start iteration.
    for k=1:num_theta
        disp(k);
        % The current projection.
        f_proj = f_projections(:, k);
        
        % The current orientation for the projection.
        estimated_orientation = ...
            Orientation(first_estimate_theta(k), first_estimate_shifts(1));
        
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

                    % Take the back projection of the two quantities. 
                    S = sliceMatrix(given_orientation, projection_parameters);
                    modified_f_proj = f_proj./(noise_estimate.^2);
                    noise_proj = ...
                        ones(projection_length, 1)./(noise_estimate.^2);

                    % Back project the projections with probability as weights. 
                    back_projection = S'*modified_f_proj;
                    disp(n/d);
                    back_projection = (n/d)*back_projection*resolution_angle*resolution_space; 

                    noise_back_projection = S'*noise_proj;
                    noise_back_projection = (n/d)*noise_back_projection*resolution_angle*resolution_space;

                    % Add the weightaged images.
                    numerator_image = numerator_image + back_projection;
                    denominator_image = denominator_image + noise_back_projection;

                    % Calculate the noise estimates.
                    c_proj = S*image_estimate;
                    for r=1:projection_length
                        next_noise_vector(r, k) = next_noise_vector(r, k) +...
                            0.5*(n/d)*(norm(c_proj(r) - f_proj(r))^2)*resolution_angle*resolution_space;
                    end
                end
            end
        end
    end

    % Calculate the next image estimate.
    next_image_estimate = ...
        numerator_image./(denominator_image + 1./mean(prior_variance(:)));

    % Reshape the fourier transform to show the image.
    img = reshape(next_image_estimate, [output_size, output_size]);
    next_image = real(ifft2(img));
    figure; imshow(next_image, []);

    % Calculate the scaled next image estimate.
    next_image_estimate = fft2(next_image);
    next_image_estimate = next_image_estimate(:);

    next_noise_vector = sqrt(next_noise_vector);

    % Estimates for the nest iteration.
    noise_estimate = mean(next_noise_vector(:));
    image_estimate = next_image_estimate;

    % Calculate the variance of the prior gaussian distribution.
    prior_variance = zeros(size(image_estimate));
    for i=1:size(image_estimate)
        prior_variance(i) = 0.5*norm(image_estimate(i))^2;
    end
    prior_variance = prior_variance*4;
end