classdef Orientation
   properties
      theta
      shift
   end
   methods
      function obj = Orientation(theta, shift)
          obj.theta = theta;
          obj.shift = shift;
      end
   end
end