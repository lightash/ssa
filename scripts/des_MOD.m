function [des, q] = des_MOD(S1, S2, X)
%des_MOD    Decision making algoritm, based on mutually 
%    orthogonal decomposition.
%    S1 and S2 are standards to which the decision is made;
%    X is input vector;
%    des is number of chosen standard (1 or 2);
%    q is value in favor of each competing version.

ls1 = len(S1);    % length of S1
ls2 = len(S2);    %           S2

es1 = nrm(S1);    % orth of S1
es2 = nrm(S2);    %         S2

ls12 = S1*es2';   % length of projection of S1 on es2
ls21 = S2*es1';   %                         S2    es1

% s12 = ls12*es2;   % projection of S1 on es2
s21 = ls21*es1;   %               S2    es1

% s11 = S1 - s12;   % exclusive component of S1
% s22 = S2 - s21;   %                        S2
% 
% ls11 = len(s11);  % length of s11
% ls22 = len(s22);  %           s22

lx = len(X);

px21 = X * es1';  % length of projection
px12 = X * es2';

x21 = px21*es1;   % projection
x12 = px12*es2;

x22 = X - x21;    % exclusive component
x11 = X - x12;

mu0 = len(s21/ls1);
mu1 = len(x11/lx);
mu2 = len(x22/lx);

dmu1 = abs(mu1 - mu0);
dmu2 = abs(mu2 - mu0);

w11 = (1 - dmu1/mu0)/3;
w22 = (1 - dmu2/mu0)/3;

x1 = x21;
x2 = x12;

mdmu1 = abs(len(x1) - ls1);
mdmu2 = abs(len(x1) - ls21);

u1 = (1 - mdmu1/abs(ls1-ls21))/3;
u2 = (1 - mdmu2/abs(ls1-ls21))/3;

mdmu1 = abs(len(x2) - ls12);
mdmu2 = abs(len(x2) - ls2);

v1 = (1 - mdmu1/abs(ls2 - ls12))/3;
v2 = (1 - mdmu2/abs(ls2 - ls12))/3;

q(1) = w11 + u1 + v1;   % decision for S1
q(2) = w22 + u2 + v2;   %              S2

if q(1) >= q(2)
   des = 1;
else
   des = 2;
end