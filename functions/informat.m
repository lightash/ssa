function prob = informat(a, b)
%informat   Calculates probabilities of
% zero and alternative hypothesis and
% errors of the first and second kind
%   prob - vector with probabilities
%   a & b - data vectors

ord = -1:.1:1;

H1 = hist(a,ord);
H2 = hist(b,ord);

H1 = H1/sum(H1);
H2 = H1/sum(H2);

[~,thr] = max(H1 + H2);

P1 = H1 ./ (H1 + H2);
P2 = H2 ./ (H1 + H2);

for i = 1:length(P1)
   if isnan(P1(i))
      if i < thr
         P1(i) = 0;
      else
         P1(i) = 1;
      end
   end
   if isnan(P2(i))
      if i < thr
         P2(i) = 0;
      else
         P2(i) = 1;
      end
   end
end

% First digit - correct variant
% Second digit - made decision
prob(1) = sum(P2(1:thr));     % 00
prob(2) = sum(P2(thr:end));   % 01
prob(3) = sum(P1(1:thr));     % 10
prob(4) = sum(P1(thr:end));   % 11

% Norming ( 00+01=1  10+11=1 )
s = prob(1) + prob(2);
prob(1:2) = prob(1:2) / s;
s = prob(3) + prob(4);
prob(3:4) = prob(3:4) / s;
