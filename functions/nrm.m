function [output, scale] = nrm(input, without_mean)
% Norming by Euclid distance.
%    without_mean is bool, whether to substrate mean value
%    from input before norming.
%    scale is a scale coefficient, length of input vector

if nargin < 2, without_mean = 0; end

output = zeros(size(input));

if ~without_mean
   for i = 1:size(input,1)
      scale = sqrt( input(i,:) * input(i,:)' );
      output(i,:) = input(i,:) ./ scale;
   end
else
   for i = 1:size(input,1)
      input(i,:) = input(i,:) - mean(input(i,:));
      scale = sqrt( input(i,:) * input(i,:)' );
      output(i,:) = input(i,:) ./ scale;
   end
end
