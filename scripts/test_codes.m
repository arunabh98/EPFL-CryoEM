% Get the image.
P = phantom(200);

% Pad the image with a fixed boundary of 3 pixels.
% P = padarray(P, [3th, 3], 0.0);

% Constants.
sigmaNoiseFraction = 0.00;
max_shift_amplitude = 0;
filename = ...
    '../results/bayesian_estimation/error_angles_and_shifts/5_percent_noise/';
num_theta = 180;
max_angle_err = 1;
max_shift_err = 0;
resolution_angle = 1;
resolution_space = 1;
no_of_iterations = 1;
mask=ones(size(P));
n = size(P, 1);
L_pad = 3032; 

% Things to write in the observation file.
theta_to_write = zeros(10, num_theta);

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

fourier_radial  = ...
    backproject_fourier_alternate(f_projections, first_estimate_theta);
first_estimate_model_alternate = Ifft2_2_Img(fourier_radial, L_pad);
figure; imshow(first_estimate_model_alternate);

modified_f_projections = zeros(size(f_projections));

error = 0;
theta(1) = theta(1);
for i=1:num_theta
%     c_proj = project_fourier_alternate(fourier_radial,...
%         theta(i), 287);
    c_proj_1 = project_fourier_alternate(fourier_radial,...
        first_estimate_theta(i), 287);
%     disp(size(c_proj, 1));
%     c_proj(isnan(c_proj)) = complex(0, 0);
%     disp(c_proj(ceil(size(c_proj, 1)/2) + 1));
    f_proj = f_projections(:, i);
%     disp(norm(c_proj - f_proj));
%     disp(norm(c_proj_1 - f_proj));
    error = error + norm(c_proj_1 - f_proj);
    modified_f_projections(:, i) = c_proj_1;
end

fourier_radial_modified = ...
    backproject_fourier_alternate(modified_f_projections, first_estimate_theta);
first_estimate_model = Ifft2_2_Img(fourier_radial_modified, L_pad);
figure; imshow(first_estimate_model);

error/num_theta

