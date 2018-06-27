% Get the image.
P = phantom(200);

% Pad the image with a fixed boundary of 3 pixels.
P = padarray(P, [3, 3], 0.0);

% Constants.
max_shift_amplitude = 4;

% Define ground truth angles and take the tomographic projection.
theta = 0:2:178;
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
theta_to_write(6, :) = original_shifts;

% Transform all entities to the frequency space.
f_projections = ifftshift(projections,1);
f_projections = fft(f_projections, [ ], 1);
f_projections = fftshift(f_projections, 1);

% Change original projections to frequency domain.
of_projections = ifftshift(original_projections,1);
of_projections = fft(of_projections, [ ], 1);
of_projections = fftshift(of_projections, 1);

first_part = (1/(295*2):1/295:0.5)';
second_part = -first_part(1:end-1);
frequencies = [first_part; flipud(second_part)];

plot(frequencies);

% Shift the projections.
for i=1:size(projections, 2)
    shift = original_shifts(i);
    if mod(shift, 2) == 0
        freq_compensation = exp(-1j*2*pi*frequencies*shift);
    else
        freq_compensation = -exp(-1j*2*pi*frequencies*shift);
    end
    comp_projections = of_projections(:, i).*freq_compensation;
    disp(norm(comp_projections - f_projections(:, i)));
end




