function res = sliceMatrix(...
    orientation, projection_parameters)
    res.orientation = orientation;
    res.projection_parameters = projection_parameters;
    res.adjoint = 0;
    
    % Register this variable as a sliceMatrix class
    res = class(res, 'sliceMatrix');