function res = mtimes(A, x)
if A.adjoint == 0 %A*x
     x = reshape(x, [A.projection_parameters.height, A.projection_parameters.width]);
     x = ifftshift(ifft2(fftshift(x)));
     projection = radon(x, A.orientation.theta);
     projection = circshift(projection, A.orientation.shift);
     res = fftshift(fft(ifftshift(projection)));
else %At*x
    projection = ifftshift(ifft(fftshift(x)));
    projection = circshift(projection, -A.orientation.shift);
    backprojection = iradon([projection projection],...
        [A.orientation.theta A.orientation.theta],...
        A.projection_parameters.output_size);

    % 2 D Fourier transformation
    reconstrution2DFT = fftshift(fft2(ifftshift(backprojection)));
    res = reconstrution2DFT(:);
end