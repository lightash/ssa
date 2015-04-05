function [dkl] = DKL(p, q)
% Kullback–Leibler divergence.
%    A non-symmetric measure of the difference 
%    between two probability distributions.

Dsize = size(p,2);
dkl = zeros(Dsize,2);
for i = 1:Dsize
   dkl(i,1) = sum( p(:,i) .* log( p(:,i)./q(:,i) ) );
   dkl(i,2) = sum( q(:,i) .* log( q(:,i)./p(:,i) ) );
end