clc;
close all;
clear all;
%{
% Processing raw data
ch_num = 64;
cf = fopen('d:\Dropbox\MATLAB\Диплом\Signals\EEG Motor Movement-Imagery Dataset\channels.txt');
for i = 1:ch_num
   channels{i} = fgetl(cf);
end
fclose(cf);

for patient = 1:10
   
   str = num2str(patient,'%.2d');
   load(['d:\Dropbox\MATLAB\Диплом\Signals\EEG Motor Movement-Imagery Dataset\S0' str 'R03_edfm.mat'])
   
   for i = 1:ch_num
      ch{i} = val(i,:);
   end
   ann = val(ch_num+1,:);
   Fd = 160;
   Ts = length(ch{1});
   j = 0;
   for i = 7:Ts
      if (ann(i)>5e3) && (ann(i-1)>5e3) && (ann(i-2)>5e3) && (ann(i-3)>5e3) && (ann(i-4)>5e3) && (ann(i-5)>5e3) && (ann(i-6)>5e3)
         j = j+1;
         anni(j) = i-6;
      end
   end
   j = 1;
   mar(1) = 1;
   for i = 2:length(anni)
      if anni(i)-anni(i-1)>500
         j = j+1;
         mar(j) = anni(i);
      end
   end
   mar = [mar Ts];

   T0 = 1:2:29;
   T1 = [4 6 12 14 18 24 26];
   T2 = [2 8 10 16 20 22 28];

   for i = 1:length(T0)
      m = mar( T0(i) ) : mar( T0(i)+1 );
      for j = 1:ch_num
         sil{j}(i,1:length(m)) = ch{j}(m);
      end
   end
   for i = 1:length(T1)
      m = mar( T1(i) ) : mar( T1(i)+1 );
      for j = 1:ch_num
         mov1{j}(i,1:length(m)) = ch{j}(m);
      end
      m = mar( T2(i) ) : mar( T2(i)+1 );
      for j = 1:ch_num
         mov2{j}(i,1:length(m)) = ch{j}(m);
      end
   end

   save(['d:\Dropbox\MATLAB\Диплом\Signals\EEG Motor Movement-Imagery Dataset\Processed\S0' str 'R03'],'ch','sil','mov1','mov2','Fd','Ts','mar','channels','ch_num')
end
   %%
%}
% The most different channel
diff = zeros(10,64);
for pat = 1:1
   
   str = num2str(pat,'%.2d');
   load(['d:\Dropbox\MATLAB\Диплом\Signals\EEG Motor Movement-Imagery Dataset\Processed\S0' str 'R03'])
   minl = min(size(mov1{1},2),size(mov2{1},2));

   for i = 1:ch_num
      for j = 1:7
         diff(pat,i) = diff(pat,i) + max(corr3( mov1{i}(j,1:minl) , mov2{i}(j,1:minl) ));
      end
   end
%    figure,subplot(211),plot(mov1{22}(j,1:minl)),subplot(212),plot(mov2{22}(j,1:minl))
%    diff(pat,:) = (diff(pat,:) - min(diff(pat,:)))/max(diff(pat,:) - min(diff(pat,:)));
   
   chch = [22 35 44];
   for i = 1:length(chch)
      for j = 1:7
         for k = 1:minl
            ch_sa{i}(k+minl*(j-1)) = sum(abs( mov1{chch(i)}(j,1:k) - mov2{chch(i)}(j,1:k) ));
            ch_mx{i}(k+minl*(j-1)) = max((corr3( mov1{chch(i)}(j,1:k) , mov2{chch(i)}(j,1:k) )));
         end
%          diff(i) = diff(i) + max((corr3( mov1{i}(j,1:minl) , mov2{i}(j,1:minl) )));
      end
   end
% for j = 1:7
%    hold on,stem(diff(pat,:))
% end
   figure,plot(ch_sa{1},'b')
   hold on,plot(ch_sa{2},'r')
   hold on,plot(ch_sa{3}','g')
   figure,plot(ch_mx{1}','b')
   hold on,plot(ch_mx{2}','r')
   hold on,plot(ch_mx{3}','g')
%    figure,plot(ch22(5,:)-ch35(5,:))
end
% figure,stem(mean(diff))
% figure,stem(std(diff))
%%
load('d:\Dropbox\MATLAB\Диплом\Signals\EEG Motor Movement-Imagery Dataset\S001R03_3236')
% load('~/Files/Dropbox/MATLAB/Р”РёРїР»РѕРј/Signals/EEG Motor Movement-Imagery Dataset/S001R03_3236')

T = 2:floor(sqrt(min(diff(mar))));
F = Fd./T;
t = 1:min(diff(mar));

% part = 1;
% % for part = 1:29
%    if part <= 15
%       curr_part = ch1.t0(part,:);
%       disp('t0')
%    else
%       if part <= 15+7
%          curr_part = ch1.t1(part-15,:);
%          disp('t1')
%       else
%          curr_part = ch1.t2(part-15-7,:); 
%          disp('t2')
%       end
%    end
%    disp(part)

sig1 = ch1.t0(1,:);
sig2 = ch1.t1(1,:);

% f = 0:160/size(ch1.t0,2):160-160/size(ch1.t0,2);
% figure
% for i = 1:size(ch1.t0,1)
%    sp(i,:) = abs(fft(ch1.t0(i,:)))/size(ch1.t0,2)*2;
%    plot(f(1:401),sp(i,1:401)),hold on
% end
% axis tight
f1 = 0:160/size(ch1.t1,2):160-160/size(ch1.t1,2);
figure(1)
for i = 1:size(ch1.t1,1)
   sp1(i,:) = abs(fft(ch1.t1(i,:)))/size(ch1.t1,2)*2;
   plot(f1(1:end/2),sp1(i,1:end/2)),hold on
end
axis tight
% f2 = 0:160/size(ch1.t2,2):160-160/size(ch1.t2,2);
% figure
% for i = 1:size(ch1.t2,1)
%    sp2(i,:) = abs(fft(ch1.t2(i,:)))/size(ch1.t2,2)*2;
%    plot(f2(1:end/2),sp2(i,1:end/2)),hold on
% end
% axis tight

% [Sil,amp] = perPort(sig1,14:18);
%    disp('Silence portrait is done.')
% t = 1:size(amp,2);
% figure,plot(t,amp(1,:),':r',t,amp(2,:),':g',t,min(amp(1,:),amp(2,:)),'.-b')
for j = 1:7
   [Mov,amp] = perPort(ch1.t1(j,:),3:11);
      disp('Movement portrait is done.')
   for i = 1:9
      bl(i) = len(Mov.period{i}.window{1}.per_win(1,:));
   end
   figure(1),hold on
   stem(160./[11:-1:3],bl(end:-1:1),'.r')
   
   figure(2),hold on
   plot(Mov.period{4}.window{1}.per_win(1,:))
end
% figure,plot(t,amp(1,:),':r',t,amp(2,:),':g',t,min(amp(1,:),amp(2,:)),'.-b')
%%
% close all

% pic 2
per = 19;
t = [1:per]/160*1e3;

figure('Color','w')
subplot(5,2,1),plot(t,sig1(1:per),'.-k'),axis tight
   title({'Action 1','Signal a1'},'FontName','Times New Roman','FontSize',10)
   ylabel('EEG a1, \muV')
subplot(5,2,3),plot(t,Sil.period{per-1}.periodic(1,:),'.-k'),axis tight,hold on,plot([1 t(end)],[0 0],':k')
   title('Continuous periodic component','FontName','Times New Roman','FontSize',10)
   ylabel('b1, \muV')
subplot(5,2,5),plot(t,Sil.period{per-1}.periodic(2,:),'.-k'),axis tight,hold on,plot([1 t(end)],[0 0],':k')
   title('Impulse periodic components','FontName','Times New Roman','FontSize',10)
subplot(5,2,7),plot(t,Sil.period{per-1}.periodic(3,:),'.-k'),axis tight,hold on,plot([1 t(end)],[0 0],':k')
subplot(5,2,9),plot(t,Sil.period{per-1}.periodic(per,:),'.-k'),axis tight
   xlabel({'Time t, ms' 'a'},'FontName','Times New Roman','FontSize',10)
   title('...................','FontSize',24)

subplot(5,2,2),plot(t,sig2(1:per),'.-k'),axis tight
   title({'Action 2','Signal a1'},'FontName','Times New Roman','FontSize',10)
   ylabel('EEG a1, \muV')
subplot(5,2,4),plot(t,Mov.period{per-1}.periodic(1,:),'.-k'),axis tight,hold on,plot([1 t(end)],[0 0],':k')
   title('Continuous periodic component','FontName','Times New Roman','FontSize',10)
   ylabel('b1, \muV')
subplot(5,2,6),plot(t,Mov.period{per-1}.periodic(2,:),'.-k'),axis tight,hold on,plot([1 t(end)],[0 0],':k')
   title('Impulse periodic components','FontName','Times New Roman','FontSize',10)
subplot(5,2,8),plot(t,Mov.period{per-1}.periodic(3,:),'.-k'),axis tight,hold on,plot([1 t(end)],[0 0],':k')
subplot(5,2,10),plot(t,Mov.period{per-1}.periodic(per,:),'.-k'),axis tight
   title('...................','FontSize',24)
   xlabel({'Time t, ms' 'b'},'FontName','Times New Roman','FontSize',10)

% pic 3
per = 24;
t = [1:per]/160*1e3;
for imp = 1:per
   int_imp1(imp) = len(Sil.period{per-1}.periodic(per-imp+1,:));
   int_imp2(imp) = len(Mov.period{per-1}.periodic(per-imp+1,:));
end
for win = 1:length(Sil.period{per-1}.window)
   int_win1(win) = len(Sil.period{per-1}.window{win}.per_win(per-19+1,:));
   int_win2(win) = len(Mov.period{per-1}.window{win}.per_win(per-19+1,:));
end

figure('Color','w')
subplot(2,2,1),stem(t,int_imp1,'.k'),axis tight
   xlabel('Length of impulses \tau_i','FontName','Times New Roman','FontSize',10)
   ylabel('Components intensity A, \muV','FontName','Times New Roman','FontSize',10)
subplot(2,2,3),hist(int_win1,20),axis tight
   xlabel('Components intensity A_1_9, \muV','FontName','Times New Roman','FontSize',10)
   ylabel('Relative frequency','FontName','Times New Roman','FontSize',10)

subplot(2,2,2),stem(t,int_imp2,'.k'),axis tight
   xlabel('Length of impulses \tau_i','FontName','Times New Roman','FontSize',10)
   ylabel('Components intensity A, \muV','FontName','Times New Roman','FontSize',10)
subplot(2,2,4),hist(int_win2,20),axis tight
   xlabel('Components intensity A_1_9, \muV','FontName','Times New Roman','FontSize',10)
   ylabel('Relative frequency','FontName','Times New Roman','FontSize',10)
%%
close all
% pic 4
t = [1:T(end)]/Fd;
f = 1./t;
int_per1(1) = abs(mean(sig1));
int_per2(1) = abs(mean(sig2));
for per = 1:T(end)-1
   int_per1(per+1) = len(Sil.period{per}.periodic(1,:));
   int_per2(per+1) = len(Mov.period{per}.periodic(1,:));
end


per = 25;
comp = per - 10;
et1 = Sil.period{per-1}.periodic(1+comp,:);
et1 = et1/len(et1);
et1 = transform(et1,'matrix_repeat');
et1(per-comp+1:end,:) = zeros(comp,per);

et2 = Mov.period{per-1}.periodic(1+comp,:);
et2 = et2/len(et2);
et2 = transform(et2,'matrix_repeat');
et2(per-comp+1:end,:) = zeros(comp,per);

etR = ones(1,per) / sqrt(per);

noise = linspace(0,2,20);
stn = 1e4;  % Statistics size
ord = -1:.1:1;  % Histogram bars

for noi = 1:length(noise)
   disp([noi length(noise)])
   
   for st = 1:stn
      
      sig = noise(noi)*randn(per) + transform( sig1(1:per^2) ,'matrix');
      
      for i = 1:per
         sig(i,:) = sig(i,:)/len(sig(i,:));
         
         R1(st,i) = sig(i,:) * et1(i,:)';
         R2(st,i) = sig(i,:) * et2(i,:)';
%          K(1,st) = sig(i,:) * et1';
%          K(2,st) = sig(i,:) * et2';
      end
      
      R1(st,:) = R1(st,:)/len(R1(st,:));
      R2(st,:) = R2(st,:)/len(R2(st,:));
      
      K(1,st) = R1(st,:) * etR';
      K(2,st) = R2(st,:) * etR';
      
   end
   
   H1(noi,:) = hist(K(1,:),ord);
   H1(noi,:) = H1(noi,:)/sum(H1(noi,:));
   H2(noi,:) = hist(K(2,:),ord);
   H2(noi,:) = H2(noi,:)/sum(H2(noi,:));
   figure,plot(H1(noi,:),'g'),hold on,plot(H2(noi,:))
   
   [~,thr(noi)] = max(H1(noi,:)+H2(noi,:));
   
   for i = 1:length(ord)
      
      P1(noi,i) = H1(noi,i)./(H1(noi,i)+H2(noi,i));
      if isnan(P1(noi,i))
         if i < thr(noi)
            P1(noi,i) = 0;
         else
            P1(noi,i) = 1;
         end
      end
      
%       P2(noi,i) = H(2,i)./(H(1,i)+H(2,i));
%       if isnan(P2(noi,i))
%          if i < thr(noi)
%             P2(noi,i) = 1;
%          else
%             P2(noi,i) = 0;
%          end
%       end
   
   end
   
   p11(noi) = sum( P1( noi, thr(noi):end ) );
   p10(noi) = sum( P1( noi, 1:thr(noi) ) );
   p11(noi) = p11(noi)/(p11(noi)+p10(noi));
   p10(noi) = p10(noi)/(p11(noi)+p10(noi));
   
%    p11(noi) = p11(noi)/sum(p11(noi));
%    p10(noi) = p10(noi)/sum(p10(noi));
   
   
end

%%
% close all
for noi = 1:10:length(noise)
%    figure,plot(H1(noi,:),'.-'),hold on,plot(H2(noi,:),'.:r')
%    figure,plot(P1(noi,:),'.-'),axis tight
end

maxf = 40;
figure('Color','w')
subplot(2,2,1),stem(f(f<maxf),int_per1(f<maxf),'.k'),axis tight
   xlabel('Frequency f, Hz','FontName','Times New Roman','FontSize',12)
   ylabel('Amplitudes A  _f ,  \muV','FontName','Times New Roman','FontSize',12)
subplot(2,2,2),stem(f(f<maxf),int_per2(f<maxf),'.k'),axis tight
   xlabel('Frequency f, Hz','FontName','Times New Roman','FontSize',12)
   ylabel('Amplitudes A  _f ,  \muV','FontName','Times New Roman','FontSize',12)
subplot(2,2,3:4),plot(p11,'k'),hold on,plot(p10,'--k'),axis tight
   xlabel('Noise intensity \sigma','FontName','Times New Roman','FontSize',12)
   ylabel('Probabilities estimates','FontName','Times New Roman','FontSize',12)
   hleg = legend('$\hat{p}$11','$\hat{p}$10');
   set(hleg,'Interpreter','latex','FontName','Times New Roman','FontSize',14)



% % Rhytms comparison
% figure('Color','w')
% subplot(3,1,1:2),plot(sig1(1:t(end)),'k')
% 
% for per = 1:length(T)
% %    for imp = 1:T(per)
%       hold on,plot(Sil.period{per-1}.per_long(1,:),'.-b')
% %    end
% end
% % hold on,plot(Port.per_long_sum,'.-r'),axis tight
% 
% for per = 1:length(T)
%    spec(per) = len(Sil.period{per}.periodic(1,:));
% end
% subplot(3,1,3),stem(T,spec,'b')
% % hold on,stem(0,len(Port.per_long_sum),'r')
% axis tight,ylabel('')

% delta = [0 4];
% teta = [4 7];
% alpha = [7 14];
% beta = [15 maxf];
% gamma = [maxf 100];
% mu = [8 13];
% 
% rhD = zeros(1,t(end));
% rhT = zeros(1,t(end));
% rhA = zeros(1,t(end));
% rhB = zeros(1,t(end));
% rhG = zeros(1,t(end));
% rhM = zeros(1,t(end));
% for i = 1:length(T)
%    if F(i)>delta(1) && F(i)<delta(2)
%       rhD(1:length(c1{i}.b)) = rhD(1:length(c1{i}.b)) + c1{i}.b;
%    end
%    if F(i)>teta(1) && F(i)<teta(2)
%       rhT(1:length(c1{i}.b)) = rhT(1:length(c1{i}.b)) + c1{i}.b;
%    end
%    if F(i)>alpha(1) && F(i)<alpha(2)
%       rhA(1:length(c1{i}.b)) = rhA(1:length(c1{i}.b)) + c1{i}.b;
%    end
%    if F(i)>beta(1) && F(i)<beta(2)
%       rhB(1:length(c1{i}.b)) = rhB(1:length(c1{i}.b)) + c1{i}.b;
%    end
%    if F(i)>gamma(1) && F(i)<gamma(2)
%       rhG(1:length(c1{i}.b)) = rhG(1:length(c1{i}.b)) + c1{i}.b;
%    end
%    if F(i)>mu(1) && F(i)<mu(2)
%       rhM(1:length(c1{i}.b)) = rhM(1:length(c1{i}.b)) + c1{i}.b;
%    end
%    
% end
% 
% figure,plot(t,ch1.sig(1:t(end)),'k',t,rhB,'m',t,rhG,'c',t,rhA,'r',t,rhT,'g',t,rhM,'b'),axis tight
% legend('F3','beta','gamma','alpha','teta','mu')
% 
% ch1.sp = fft(ch1.sig(1:t(end)));
% figure,stem(abs(ch1.sp(1:end)),'.k'),axis tight
% 
% for i = 1:length(T)
%    c1{i}.sp = fft(c1{i}.b);
%    hold on,stem(abs(c1{i}.sp(1:end)),'.b')
% end















