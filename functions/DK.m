function [J] = DK(p, q, N1, N2)
% Kullback difference criterium.
%    p and q are hist_len by num_hists matrices, built on N1 and N2
%    observations accordingly.
%    J is 1 by num_hists vector of symmetric Kullback difference numbers.

num_hists = size(p,2);

NN = N1*N2/(N1+N2);

J = zeros(1,num_hists);
for i = 1:num_hists
   J(1,i) = NN*sum( (p(:,i)-q(:,i)) .* log( p(:,i)./q(:,i) ) );
end