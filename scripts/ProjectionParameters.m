classdef ProjectionParameters
   properties
      width
      height
      output_size
      projection_length
   end
   methods
      function obj = ...
              ProjectionParameters(width, height, output_size, projection_length)
          obj.width = width;
          obj.height = height;
          obj.output_size = output_size;
          obj.projection_length = projection_length;
      end
   end
end