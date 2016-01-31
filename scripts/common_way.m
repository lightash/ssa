clc;
% close all;
clear all;

load('D:\Dropbox\Signals\incartdb\I20\I20proc.mat')
in = val(1,:);
annot(2461) = 'N';   % Fusion of ventricular and normal beat

all_beats = 1:length(annot);
% Normal beat, Atrial premature beat, Premature ventricular contraction
bmark = 'NAV';
btypeN = 3;       % Beat types to examine
for i = 1:btypeN
   Bnum{i} = all_beats(annot == bmark(i));                           % anniNAV
   Blen(i) = length(Bnum{i});
   Bpos{i} = mark(Bnum{i});                                          % ann
   Bord(all_beats(annot ==  bmark(i))) = i*ones(1,length(Bpos{i}));  % annNAV
end
Bwin = [-47 80];   % Borders of PQRST period
winL = Bwin(2)-Bwin(1)+1;

% Generating portraits
disp('Generating portraits')
perN = all_beats(end);  % Number of periods to use

f = zeros(perN, winL );
in_stab = zeros(size(in));
for per = 1:perN

   period = mark(per);
   window = period+Bwin(1): period+Bwin(2);
   in_stab(window) = in(window) - mean(in(window));
   
   f(per,:) = in(window);
   f(per,:) = f(per,:) - mean(f(per,:));
%    f(per,:) = nrm(f(per,:),1);
end

port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = mean(f(Bnum{btype},:));
   port{btype} = nrm(port{btype},1);
end
[~,~,~,d] = des_MOD(port{1},port{2},f(per,:));

x = (Bwin(1):Bwin(2))/Fd*1e3;
% figure('Color','w')
% plot(x,port{2},'k'),hold on
% plot(x,d.s22,'-.k','LineWidth',1.5),hold on
% legend('S_2','s_2_2'),grid on,axis tight
% xlabel('¬рем€ от R-зубца, мс','FontName','Times New Roman','FontSize',12)
% ylabel('¬заимноортогональные составл€ющие','FontName','Times New Roman','FontSize',12)
% %%
% r = 1287;%randi(perN);
% [des,~,~,d] = des_MOD(port{1},port{2},f(r,:));
% figure('Color','w')
% subplot(121),stem([d.w11 d.u1 d.v1],'.k'),axis([0 4 -.3 .3]),grid
% set(gca,'XTickLabel','')
% title(Bord(r))
% ylabel('«начение критери€')
% % ylabel(r)
% xlabel('w11      u1      v1')
% subplot(122),stem([d.w22 d.u2 d.v2],'.k'),axis([0 4 -.3 .3]),grid
% set(gca,'XTickLabel','')
% xlabel('w22      u2      v2')
% % title(des)

% figure('Color','w')
% plot(x,port{1},'k',x,port{2},'--k',x,d.s12,'-.k',x,d.s21,':k'),hold on
% plot(x,d.s11,'-.k',x,d.s22,':k','LineWidth',1.5),hold on
% legend('S_1','S_2','s_1_2','s_2_1','s_1_1','s_2_2'),grid on,axis tight
% % xlabel('¬рем€ от R-зубца, мс','FontName','Times New Roman','FontSize',12)
% % ylabel('¬заимноортогональные составл€ющие','FontName','Times New Roman','FontSize',12)
% xlabel('„ас в≥дносно R-зубц€, мс','FontName','Times New Roman','FontSize',12)
% ylabel('¬заЇмно-ортогональн≥ складов≥','FontName','Times New Roman','FontSize',12)
% figure('Color','w')
% plot(x,port{1},'k',x,port{2},'--k'),hold on
% legend('Normal beats','Atrial premature beats'),grid on,axis tight
% xlabel('¬рем€ от R-зубца, мс','FontName','Times New Roman','FontSize',12)
% ylabel('‘орма QRS-комплексов по типам','FontName','Times New Roman','FontSize',12)
% xlabel('„ас в≥дносно R-зубц€, мс','FontName','Times New Roman','FontSize',12)
% ylabel('‘орма QRS-комплекс≥в по типам','FontName','Times New Roman','FontSize',12)

%%
% Guessing
disp('Guessing')
btypeN = 2;

des = zeros(btypeN);
for per = 1:perN
   disp(per)
   for btype = 1:btypeN
      cor(btype,per) = f(per,:) * port{btype}';
      cor(btype,per) = (cor(btype,per) +1)/2;
   end
   [~,ind] = max(cor(:,per));
      
      [~,q(per,:)] = des_MOD(port{1},port{2},f(per,:));
%       [~,ind] = max(q(per,:));
      
   if ~any(per == Bnum{3})
      des(Bord(per),ind) = des(Bord(per),ind) + 1/Blen(Bord(per));
   end
end

figure
k = 0;
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
      subplot(btypeN,btypeN,k),stem(des(i,j),'.-'),axis([0 2 0 1])
      xlabel(des(i,j))
  end
end
% title((des(1,1)+des(2,2)+des(3,3))/3)
title((des(1,1)+des(2,2))/2)






