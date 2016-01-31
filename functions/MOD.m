function [dec] = MOD(S1, S2, X)
%MOD    Mutually orthogonal decomposition.
%    S1 and S2 are standards to which the decision is made;
%    X is input vector;
%    dec is a structure with the result of 
%    decomposition of S1, S2 and X.

dec.S1 = S1;
dec.S2 = S2;

dec.es1 = nrm(S1);         % orth of S1
dec.es2 = nrm(S2);         %         S2

ls12 = S1 * dec.es2';      % length of projection of S1 on es2
ls21 = S2 * dec.es1';      %                         S2    es1

dec.s12 = ls12 * dec.es2;  % projection of S1 on es2
dec.s21 = ls21 * dec.es1;  %               S2    es1

dec.s11 = S1 - dec.s12;    % exclusive component of S1
dec.s22 = S2 - dec.s21;    %                        S2

if nargin > 2
   
   dec.X = X;

   dec.x1 = (X * dec.es1') * dec.es1;  % projection of X on S1
   dec.x2 = (X * dec.es2') * dec.es2;  %                    S2

   % x22 = X - x1;
   % x11 = X - x2;
   dec.x22 = (X * nrm(dec.s22)') * nrm(dec.s22);   % exclusive components of X
   dec.x11 = (X * nrm(dec.s11)') * nrm(dec.s11);

   y = dec.x1 - dec.x2;
   yx11 = (dec.x11 * nrm(y)') * nrm(y);
   dec.z = len(dec.x11 - yx11);  % distance from S1-O-S2 hyperplane to X
   
end