function [I, Ia] = ED(p, q, mode)
% Energy and Dispersion informativity criterium.
%    p and q are (num,len) matrices of data.
%    mode is number of used formulae, where '*a' is
%    alternative decision:
%    (1) I = E(i) - D(i);
%    (2) I = E(i) - D(i) - Ea(i);
%    (3) I = E(i) - D(i) - Ea(i) - Da(i).
%    I and Ia are (1,len) vectors of informativity numbers.

if nargin == 2, mode = 3; end
if nargin == 1, mode = 1; end

len_p = size(p,2);

I = zeros(1,len_p);
E = zeros(1,len_p);
D = zeros(1,len_p);
if nargin > 1
   Ia = zeros(1,len_p);
   Ea = zeros(1,len_p);
   Da = zeros(1,len_p);
end

for i = 1:len_p
   E(i) = mean(p(:,i))^2;
   D(i) = std(p(:,i))^2;
   I(i) = E(i) - D(i);
end

if nargin > 1
   for i = 1:len_p
      Ea(i) = mean(q(:,i))^2;
      Da(i) = std(q(:,i))^2;
      Ia(i) = Ea(i) - Da(i);
      if mode > 1
         I(i) = I(i) - Ea(i);
         Ia(i) = Ia(i) - E(i);
         if mode > 2
            Da(i) = std(q(:,i))^2;
            I(i) = I(i) - Da(i);
            Ia(i) = Ia(i) - D(i);
         end
      end
   end
end