function [des, q, dec] = des_MOD(dec, win)
%des_MOD    Decision making algoritm, based on mutually 
%    orthogonal decomposition data.
%    dec should contain S1 and S2 standards and their decompositions,
%    to which the decision is made;
%    win{1:3} is a structure with a lists of counts to use for each
%    of components;
%    des is number of chosen standard (1 or 2);
%    q is value in favor of each competing version;
%    dec is a same structure to which decision data is added.

if nargin > 1

   ls1 = len(dec.S1(win{1}));
   ls2 = len(dec.S2(win{2}));
   
   ls21 = dec.S2(win{1}) * dec.es1(win{1})';
   ls12 = dec.S1(win{2}) * dec.es2(win{2})';

   % u - collinear 1
   dmu1_u = abs(len(dec.x1(win{1})) - ls1);
   dmu2_u = abs(len(dec.x1(win{1})) - ls21);

   dec.u1 = (1 - dmu1_u/abs(ls1-ls21))/3;
   dec.u2 = (1 - dmu2_u/abs(ls1-ls21))/3;

   % v - collinear 2
   dmu1_v = abs(len(dec.x2(win{2})) - ls12);
   dmu2_v = abs(len(dec.x2(win{2})) - ls2);

   dec.v1 = (1 - dmu1_v/abs(ls2 - ls12))/3;
   dec.v2 = (1 - dmu2_v/abs(ls2 - ls12))/3;

   % w - exclusive
   mu0 = len(dec.s21(win{3})/len(dec.S1(win{3})));
   mu1 = len(dec.x11(win{3})/len(dec.X(win{3})));
   mu2 = len(dec.x22(win{3})/len(dec.X(win{3})));

   dmu1_w = abs(mu1 - mu0);
   dmu2_w = abs(mu2 - mu0);

   dec.w11 = (1 - dmu1_w/mu0)/3;
   dec.w22 = (1 - dmu2_w/mu0)/3;

else  % without window data

   ls1 = len(dec.S1);  % length of S1
   ls2 = len(dec.S2);  %           S2

   ls21 = dec.S2 * dec.es1';  % length of projection of S2 on es1
   ls12 = dec.S1 * dec.es2';  %                         S1    es2

   % u - collinear 1
   dmu1_u = abs(len(dec.x1) - ls1);
   dmu2_u = abs(len(dec.x1) - ls21);

   dec.u1 = (1 - dmu1_u/abs(ls1-ls21))/3;
   dec.u2 = (1 - dmu2_u/abs(ls1-ls21))/3;

   % v - collinear 2
   dmu1_v = abs(len(dec.x2) - ls12);
   dmu2_v = abs(len(dec.x2) - ls2);

   dec.v1 = (1 - dmu1_v/abs(ls2 - ls12))/3;
   dec.v2 = (1 - dmu2_v/abs(ls2 - ls12))/3;

   % w - exclusive
   mu0 = len(dec.s21/len(dec.S1));
   mu1 = len(dec.x11/len(dec.X));
   mu2 = len(dec.x22/len(dec.X));

   dmu1_w = abs(mu1 - mu0);
   dmu2_w = abs(mu2 - mu0);

   dec.w11 = (1 - dmu1_w/mu0)/3;
   dec.w22 = (1 - dmu2_w/mu0)/3;

end

q(1) = dec.w11 + dec.u1 + dec.v1;   % decision for S1
q(2) = dec.w22 + dec.u2 + dec.v2;   %              S2

[~,des] = max(q);
