function [des, q, z, d] = des_MOD(S1, S2, X)
%des_MOD    Decision making algoritm, based on mutually 
%    orthogonal decomposition.
%    S1 and S2 are standards to which the decision is made;
%    X is input vector;
%    des is number of chosen standard (1 or 2);
%    q is value in favor of each competing version;
%    z is a distance from S1-O-S2 hyperplane to X;
%    d is a structure with decomposed vectors of S1, S2, X.

ls1 = len(S1);          % length of S1
ls2 = len(S2);          %           S2

d.es1 = nrm(S1);        % orth of S1
d.es2 = nrm(S2);        %         S2

ls12 = S1 * d.es2';     % length of projection of S1 on es2
ls21 = S2 * d.es1';     %                         S2    es1

d.s12 = ls12 * d.es2;   % projection of S1 on es2
d.s21 = ls21 * d.es1;   %               S2    es1

d.s11 = S1 - d.s12;     % exclusive component of S1
d.s22 = S2 - d.s21;     %                        S2

d.x1 = (X * d.es1') * d.es1;  % projection of X on S1
d.x2 = (X * d.es2') * d.es2;  %                    S2

% x22 = X - x1;
% x11 = X - x2;
d.x22 = (X * nrm(d.s22)') * nrm(d.s22);   % exclusive components of X
d.x11 = (X * nrm(d.s11)') * nrm(d.s11);

mu0 = len(d.s21/ls1);
mu1 = len(d.x11/len(X));
mu2 = len(d.x22/len(X));

dmu1_w = abs(mu1 - mu0);
dmu2_w = abs(mu2 - mu0);

d.w11 = (1 - dmu1_w/mu0)/3;
d.w22 = (1 - dmu2_w/mu0)/3;

dmu1_u = abs(len(d.x1) - ls1);
dmu2_u = abs(len(d.x1) - ls21);

d.u1 = (1 - dmu1_u/abs(ls1-ls21))/3;
d.u2 = (1 - dmu2_u/abs(ls1-ls21))/3;

dmu1_v = abs(len(d.x2) - ls12);
dmu2_v = abs(len(d.x2) - ls2);

d.v1 = (1 - dmu1_v/abs(ls2 - ls12))/3;
d.v2 = (1 - dmu2_v/abs(ls2 - ls12))/3;

q(1) = d.w11 + d.u1 + d.v1;   % decision for S1
q(2) = d.w22 + d.u2 + d.v2;   %              S2

if q(1) >= q(2)
   des = 1;
else
   des = 2;
end

y = d.x1 - d.x2;
yx11 = (d.x11 * nrm(y)') * nrm(y);
z = len(d.x11 - yx11);