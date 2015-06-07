function [normalized] = nrm(input, without_mean)
% Norming by Euclid distance.
%    without_mean is bool, whether to substrate mean value
%    from input before norming.

if nargin < 2, without_mean = 0; end

if ~without_mean
   normalized = input ./ sqrt( input * input' );
else
   input = input - mean(input);
   normalized = input ./ sqrt( input * input' );
end
