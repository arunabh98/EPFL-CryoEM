% Get the image.
P = phantom(200);

% Reduce the size of the image for speed.
P = imresize(P, 0.2);

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.05;
max_shift_amplitude = 1;
num_theta = 180;
max_angle_err = 0;
max_shift_err = 1;
resolution_angle = 1;
resolution_space = 1;
mask=ones(size(P));
n = size(P, 1);
L_pad = 288; 

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

% Define ground truth angles and take the tomographic projection.
% theta = datasample(0:0.5:359.5, num_theta);  
theta = 0:1:179;
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

% Initialize parameters needed for searching in the space.
prior_parameters = PriorParameters(max_angle_err, max_shift_err,...
    resolution_angle, resolution_space);

% Add noise to projections.
[projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Transform all entities to the frequency space.
f_p = ifftshift(projections,1);
f_projections = fft(f_p,[ ],1);
f_projections = fftshift(f_projections,1); % put DC central after filtering 

% In the first case the angles and shifts will be unknown upto a 
% certain limit.
first_estimate_theta = mod(theta +...
    randi([-max_angle_err, max_angle_err], 1, num_theta), 180);
first_estimate_shifts = original_shifts +...
    randi([-max_shift_err, max_shift_err], 1, num_theta);

% Begin estimation of the first model.
max_angle_err = 0;
max_shift_err = 5;
prob_matrix_height = (2*max_angle_err)/resolution_angle + 1;
prob_matrix_width = 2*max_shift_err/resolution_space + 1;
prob_matrix = ...
    zeros(prob_matrix_height, prob_matrix_width,...
        size(f_projections, 2)) + 1/(prob_matrix_height*prob_matrix_width);

% Start estimating the image.
fourier_radial = zeros(621, 621);
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

f_image_estimate = fourier_radial(:);
first_estimate_model = Ifft2_2_Img(fourier_radial, L_pad);
max_angle_err = 0;
max_shift_err = 1;

correct_shift = zeros(size(original_shifts));

error = 0;
for i=1:num_theta
    c_proj_1 = project_fourier_alternate(fourier_radial,...
        first_estimate_theta(i), first_estimate_shifts(i), 69);
    c_proj_2 = project_fourier_alternate(fourier_radial,...
        first_estimate_theta(i), first_estimate_shifts(i) + 1, 69);
    c_proj_3 = project_fourier_alternate(fourier_radial,...
        first_estimate_theta(i) - 1, first_estimate_shifts(i) - 1, 69);
    f_proj = f_projections(:, i);
    
    projections = [c_proj_3 c_proj_1 c_proj_2];
    projections_dist = bsxfun(@minus, projections, f_proj);
    [~, least_index] = min(vecnorm(projections_dist));
    correct_shift(i) = first_estimate_shifts(i)  + least_index - 2;
    
end

disp(norm(correct_shift - original_shifts, 1));
disp(norm(first_estimate_shifts - original_shifts, 1));
