function [e, c, a2] = AM(a1)
%AM    Another, easier, variant of getting periodic component from a
%   matrix.
%   [E(N), C(M), A2(M,N) ] = AM( A1(M,N) ) returns the universal periodic
%   component E and its amplitudes along axes C of a matrix A1 and remains 
%   from them - A2.

% Number of rows
M = size(a1,1);

% Centre of dots' "cloud"
b = mean(a1,1);

% Distance to centre
lb = len(b);

% Direction to the centre
e = b / lb;

% Lengths of projections on central vector
c = zeros(1,M);
for i = 1:M
   c(i) = a1(i,:) * e';
end

% Remains from substraction
a2 = zeros(size(a1));
for i = 1:M
   a2(i,:) = a1(i,:) - c(i)*e;
end

% % Leaving in 'c' only that differs from mean and giving 'e' length of that
% % mean
% mc = mean(c);
% e = e * mc;
% for i = 1:M
%    c(i) = c(i) - mc;
% end