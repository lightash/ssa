clc;
close all;
clear all;

N = 16;
x = 215:5:290;
k(1,:) = [1.47 1.08 .77 .53  .37 .32 .31 .33   .35 .37 .4  .35  .3  .26 .2 .16];
k(2,:) = [1.2  1.03 .88 .78  .72 .71  .8 .93  1.03  1  .9  .76  .68 .47 .3 .21];
n = 1:N;

S1 = k(1,:);
S2 = k(2,:);

lS1 = len(S1);
eS1 = S1/lS1;
lS2 = len(S2);
eS2 = S2/lS2;

s12 = (S1 * eS2') * eS2;
s11 = S1 - s12;
s21 = (S2 * eS1') * eS1;
s22 = S2 - s21;

mu0 = len(s11)/lS1;

stn = 2e3;  % Statistics size
ord = 0:.01:1;  % Histogram bars

noise = linspace(0,1,20);
for noi = 1:length(noise)

   ending = ['st';'nd'];
   disp(['Noise is ' num2str(noise(noi)) ' with ' num2str(1) num2str(ending(1,:)) ' as basis'])

   for st = 1:stn
      
      X = noise(noi)*randn(1,N) + S1;
      lX = len(X);
      X = X/lX;

      x21 = (X * eS1') * eS1;
      x22 = X - x21;
      x12 = (X * eS2') * eS2;
      x11 = X - x12;

      mu1 = len(x11)/lX;
      mu2 = len(x22)/lX;

      dmu1 = abs(mu0 - mu1);
      dmu2 = abs(mu0 - mu2);

      w11 = (1 - dmu1/mu0)/3;
      w22 = (1 - dmu2/mu0)/3;

      v0 = abs(len(S2)-len(s12));
      v1 = ( 1 - dmu1/v0 )/3;
      v2 = ( 1 - dmu2/v0 )/3;

      u0 = abs(len(S1)-len(s21));
      u1 = ( 1 - dmu1/u0 )/3;
      u2 = ( 1 - dmu2/u0 )/3;

      q1(st) = w11 + u1 + v1;
      q2(st) = w22 + u2 + v2;

%       disp( ['   X is S1 for ' num2str( q1*100 ) '%'] )
%       disp( ['   X is S2 for ' num2str( q2*100 ) '%'] )
      
   end
   
   H1(noi,:) = hist(q1,ord);
   H1(noi,:) = H1(noi,:)/sum(H1(noi,:));
   H2(noi,:) = hist(q2,ord);
   H2(noi,:) = H2(noi,:)/sum(H2(noi,:));
%    figure,plot(H1(noi,:)),hold on,plot(H2(noi,:),'r')
   
   [~,m1] = max(H1(noi,:));
   [~,m2] = max(H2(noi,:));
   thr(noi) = round( (m1+m2)/2 );
%    if m1>m2
%       m = m2:m1;
%    else
%       m = m1:m2;
%    end
%    thr(noi) = round(mean(find(H1(noi,m)-H2(noi,m) < 1e-1)));

%    [~,thr(noi)] = max(H1(noi,:)+H2(noi,:));

   for i = 1:length(ord)
      
      P1(noi,i) = H1(noi,i)./(H1(noi,i)+H2(noi,i));
      if isnan(P1(noi,i))
         if i < thr(noi)
            P1(noi,i) = 0;
         else
            P1(noi,i) = 1;
         end
      end
   end
   
   p11(noi) = sum( P1( noi, thr(noi):end ) );
   p10(noi) = sum( P1( noi, 1:thr(noi) ) );
   s = p11(noi)+p10(noi);
   p11(noi) = p11(noi)/s;
   p10(noi) = p10(noi)/s;
   
end
%%
figure('Color','w')
plot(noise,p11,'k'),hold on,plot(noise,p10,'--k'),axis([0 noise(end) 0 1])
xlabel('Noise intensity \sigma','FontName','Times New Roman','FontSize',12)
ylabel('Probabilities estimates','FontName','Times New Roman','FontSize',12)
hleg = legend('$\hat{p}$11','$\hat{p}$10');
set(hleg,'Interpreter','latex','FontName','Times New Roman','FontSize',14)






