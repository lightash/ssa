function [Az] = AlphaZ(p, q)
% Correlation difference criterium.
%    p and q are hist_len by num_hists matrices.
%    Az is 1 by num_hists vector of symmetric correlation difference 
%    numbers from 0 to 1.

num_hists = size(p,2);

Az = zeros(1,num_hists);
for i = 1:num_hists
   p(:,i) = nrm(p(:,i)')';
   q(:,i) = nrm(q(:,i)')';
   Az(1,i) = 1 - (p(:,i)'*q(:,i) +1)/2;
end