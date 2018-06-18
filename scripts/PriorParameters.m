classdef PriorParameters
   properties
      max_angle_err
      max_shift_err
      resolution_angle
      resolution_space
   end
   methods
      function obj = PriorParameters(...
              max_angle_err, max_shift_err, resolution_angle, resolution_space)
          obj.max_angle_err = max_angle_err;
          obj.max_shift_err = max_shift_err;
          obj.resolution_angle = resolution_angle;
          obj.resolution_space = resolution_space;
      end
   end
end