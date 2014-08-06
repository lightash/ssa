clc;
close all;
clear all;

load('d:\Dropbox\MATLAB\ƒиплом\Signals\EEG Motor Movement-Imagery Dataset\S001R03_3236')

T = 24;%2:floor(sqrt(min(diff(mar))));
F = Fd./T;
t = 1:T^2;

sig1 = ch1.t0(1,1:T^2);
sig2 = ch1.t2(1,1:T^2);

[Sil] = perPort(sig1,T);
   disp('Silence portrait is done.')

[Mov] = perPort(sig2,T);
   disp('Movement portrait is done.')

a = Sil.period{1}.window{1}.per_win(1,:);
b = Mov.period{1}.window{1}.per_win(1,:);
a = a/len(a);
b = b/len(b);
c = a;

d = abs(a-b);

noise = linspace(0,10,30);
stn = 1e3;  % Statistics size

for noi = 1:length(noise)
   disp([noi length(noise)])
   
   for st = 1:stn
      sig = noise(noi)*randn(1,T) + c;
      sig = sig/len(sig);

      R1(st) = sig * a';
      R2(st) = sig * b';
   end
   
   r(noi) = 0;
   for i = 1:length(R1)
      if R1(i)>R2(i)
         r(noi) = r(noi)+1;
      end
   end
   w(noi) = stn - r(noi);
   
   r(noi) = r(noi)/stn;
   w(noi) = w(noi)/stn;

end

tp = t/160;
tpm = [1:T]/160*1e3;

% -------------------------------
sig12 = ch1.t0(2,1:T^2);
sig22 = ch1.t2(2,1:T^2);

[Sil] = perPort(sig12,T);
   disp('Silence portrait is done.')

[Mov] = perPort(sig22,T);
   disp('Movement portrait is done.')

a2 = Sil.period{1}.window{1}.per_win(1,:);
b2 = Mov.period{1}.window{1}.per_win(1,:);
a2 = a2/len(a2);
b2 = b2/len(b2);

d2 = abs(a2-b2);

figure('Color','w')
subplot(121),plot(tp,sig1,'k'),axis tight
   ylabel('ЁЁ√, мк¬')
   xlabel({'¬рем€, с' 'а'},'FontName','Arial','FontSize',10)
subplot(122),plot(tp,sig2,'k'),axis tight
   ylabel('ЁЁ√, мк¬')
   xlabel({'¬рем€, с' 'б'},'FontName','Arial','FontSize',10)

figure('Color','w')
subplot(311),stem(tpm,a,'.k'),axis tight
   ylabel('b, мк¬')
   xlabel({'¬рем€, мс' 'а'},'FontName','Arial','FontSize',10)
subplot(312),stem(tpm,b,'.k'),axis tight
   ylabel('b, мк¬')
   xlabel({'¬рем€, мс' 'б'},'FontName','Arial','FontSize',10)
subplot(313),stem(tpm,d,'.k'),axis tight
   ylabel('–азница, мк¬')
   xlabel({'¬рем€, мс' 'в'},'FontName','Arial','FontSize',10)

% d = abs(a+a2-b-b2)/2;
%%
% ----------------------------------------------------------
% figure
% subplot(421),stem(a,'.k'),axis tight
% subplot(422),stem(b,'.k'),axis tight
% subplot(4,2,3:4),stem(d,'.k'),axis tight
% subplot(4,2,5:8),plot(noise,r,'.-k'),hold on,plot(noise,w,'.:k'),axis tight,legend('right','wrong');


[s, ix] = sort(d,'descend');

xn = ix(end);
xn1 = ix(end-1);
xn2 = ix(end-2);
xn3 = ix(end-3);
xn4 = ix(end-4);

lgc = ix(ix ~= xn & ix ~= xn1 & ix ~= xn2 & ix ~= xn3 & ix ~= xn4);
c1 = a( lgc );
a1 = a( lgc );
b1 = b( lgc );

a1 = a1/len(a1);
b1 = b1/len(b1);

for noi = 1:length(noise)
   disp([noi length(noise)])
   
   for st = 1:stn
      sig = noise(noi)*randn(1, length(c1)) + c1;
      sig = sig/len(sig);

      R1(st) = sig * a1';
      R2(st) = sig * b1';
   end
   
   r1(noi) = 0;
   for i = 1:length(R1)
      if R1(i)>R2(i)
         r1(noi) = r1(noi)+1;
      end
   end
   w1(noi) = stn - r1(noi);
   
   r1(noi) = r1(noi)/stn;
   w1(noi) = w1(noi)/stn;

end

% figure
% subplot(421),stem(a1,'.k'),axis tight
% subplot(422),stem(b1,'.k'),axis tight
% subplot(4,2,5:8),plot(noise,r1,'.-k'),hold on,plot(noise,w1,'.:k'),axis tight,legend('right','wrong');



xn = ix(1);
xn1 = ix(2);
xn2 = ix(3);
xn3 = ix(4);
xn4 = ix(5);

lgc = ix(ix ~= xn & ix ~= xn1 & ix ~= xn2 & ix ~= xn3 & ix ~= xn4);
c2 = a( lgc );
a2 = a( lgc );
b2 = b( lgc );

a2 = a2/len(a2);
b2 = b2/len(b2);

% figure
% subplot(421),stem(a2,'.k'),axis tight
% subplot(422),stem(b2,'.k'),axis tight

for noi = 1:length(noise)
   disp([noi length(noise)])
   
   for st = 1:stn
      sig = noise(noi)*randn(1, length(c2)) + c2;
      sig = sig/len(sig);

      R1(st) = sig * a2';
      R2(st) = sig * b2';
   end
   
   r2(noi) = 0;
   for i = 1:length(R1)
      if R1(i)>R2(i)
         r2(noi) = r2(noi)+1;
      end
   end
   w2(noi) = stn - r2(noi);
   
   r2(noi) = r2(noi)/stn;
   w2(noi) = w2(noi)/stn;

end

% subplot(4,2,5:8),plot(noise,r2,'.-k'),hold on,plot(noise,w2,'.:k'),axis tight,legend('right','wrong');


figure('Color','w')
plot(noise,r,'-k',noise,r1,'--k',noise,r2,':k'),hold on
plot(noise,w,'-k',noise,w1,'--k',noise,w2,':k')
   xlabel('»нтенсивность шума \sigma','FontName','Arial','FontSize',10)
   ylabel('ќценки веро€тности','FontName','Arial','FontSize',10)
   hleg = legend('ѕолный набор признаков','ќтсутствуют похожие компоненты','ќтсутствуют контрастные компоненты');










