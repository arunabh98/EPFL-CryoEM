function projection = project_fourier(fourier_radial, theta)

	% theta = theta*pi/180;
	% center = ceil(size(fourier_radial, 1)/2) + 1;
	% radius = ceil(size(fourier_radial, 1)/2)*5;

	% x1 = radius*cos(theta) + center;
	% y1 = radius*sin(theta) + center;
	% x2 = -radius*cos(theta) + center;
	% y2 = -radius*sin(theta) + center;

	% first_half = improfile(fourier_radial, [center, x1], [center, y1]);
	% first_half = first_half(1:144);
	% second_half = improfile(fourier_radial, [center, x2], [center, y2]);
	% second_half = second_half(1:145);
	% second_half = second_half(2:end);

	% second_half = flipud(second_half);

	% projection = [second_half; first_half];
	% projection(isnan(projection)) = 0;


	% radius = ceil(size(fourier_radial, 1)/2);
	% theta = theta*pi/180;
	% x1 = radius*cos(theta) + 145;
	% y1 = radius*sin(theta) + 145;
	% x2 = -radius*cos(theta) + 145;
	% y2 = -radius*sin(theta) + 145;
    image_radial = rotateAround(fourier_radial, 145, 145, theta - 90, 'bicubic');
 %    projection = improfile(fourier_radial, [x1, x2], [y1, y2], 288);
    % center = ceil(size(projection, 1)/2) + 1;
    % projection = projection(center - 144:center+143);
    % imshow(image_radial);
    projection = image_radial(:, 145);
end