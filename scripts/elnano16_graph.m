clc
close all
clear all

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
perN = all_beats(end);  % Number of periods to use
f = zeros(perN, winL);
for per = 1:perN

   period = mark(per);
   window = period+Bwin(1): period+Bwin(2);
   f(per,:) = in(window);
   fm(per) = mean(f(per,:));
   f(per,:) = f(per,:) - fm(per);

   [f(per,:), fs(per)] = nrm(f(per,:));
end
port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(Bnum{btype},:),1));
end

%% pic 1
x = (-47:80)/Fd*1e3;

figure('color','white')
plot(x,port{1},'k','LineWidth',1.5),hold on
plot(x,port{2},'--k','LineWidth',1.5),hold on
plot(x,port{3},':k','LineWidth',1.5),axis tight
leg = legend('Normal beats','Atrial premature beats','Premature ventricular contraction');
set(leg,'FontName','Times New Roman','FontSize',14,'Location','NorthOutside')
% xlabel('Time from R-peak, ms','FontName','Times New Roman','FontSize',14)
% ylabel('Shape of QRS-complexes','FontName','Times New Roman','FontSize',14),grid
xlabel('Время от R-зубца, мс','FontName','Times New Roman','FontSize',14)
ylabel('Форма QRS-комплексов','FontName','Times New Roman','FontSize',14),grid

%% pic 3
btypeN = 2;

Nbins = 16;
for i = 1:winL
   [hyp(i,:),scale(i,:)] = hist( [f(Bnum{1},i);f(Bnum{2},i)] ,Nbins);
   
   for btype = 1:btypeN
      H{btype}(:,i) = hist( f(Bnum{btype},i) ,scale(i,:))/Blen(btype)';
      % Распределить наименьший столбец по нулям
      [minh,indh] = min(H{btype}( H{btype}(:,i)>0, i));
      mzh = minh/( length( H{btype}( H{btype}(:,i)==0, i)) +1);
      H{btype}( H{btype}(:,i)==0 | H{btype}(:,i)==minh, i) = mzh;
   end
end

i = 20;
Hb = zeros(4,size(H{1}(:,i),1)); % Nback Aback Nfore Afore
for j = 1:size(H{1}(:,i),1)
   if H{1}(j,i) > H{2}(j,i)
      Hb(1,j) = H{1}(j,i);
      Hb(4,j) = H{2}(j,i);
   else
      Hb(3,j) = H{1}(j,i);
      Hb(2,j) = H{2}(j,i);
   end
end
figure('Color','w')
o = bar(scale(i,:),Hb(1,:));hold on
set(o,'FaceColor',[.5 .5 .5]);
bar(scale(i,:),Hb(2,:),'w'),hold on
o = bar(scale(i,:),Hb(3,:));hold on
set(o,'FaceColor',[.5 .5 .5]);
bar(scale(i,:),Hb(4,:),'w'),grid
legend('N','A')
% xlabel('Shape characteristic value','FontName','Times New Roman','FontSize',14)
% ylabel({'Relative frequency for' 'the intervals of values'},'FontName','Times New Roman','FontSize',14)
xlabel('Значение характеристики формы','FontName','Times New Roman','FontSize',14)
ylabel({'Относительная частота' 'интервалов значений'},'FontName','Times New Roman','FontSize',14)

%% pic 4 & 5
I(1,:) = DK( H{1}, H{2}, Blen(1), Blen(2) );
figure('Color','w')
plot(I(1,:),'.-k'),grid,axis tight
% xlabel('Sample number','FontName','Times New Roman','FontSize',14)
% ylabel('Informativity','FontName','Times New Roman','FontSize',14)
xlabel('Номер отсчёта','FontName','Times New Roman','FontSize',14)
ylabel('Информативность','FontName','Times New Roman','FontSize',14)

I(2,:) = AlphaZ( H{1}, H{2} );
figure('Color','w')
plot(I(2,:),'.-k'),grid,axis tight
% xlabel('Sample number','FontName','Times New Roman','FontSize',14)
% ylabel('Informativity','FontName','Times New Roman','FontSize',14)
xlabel('Номер отсчёта','FontName','Times New Roman','FontSize',14)
ylabel('Информативность','FontName','Times New Roman','FontSize',14)

%% pic 5 & 7
for form = 1:size(I,1)
   disp(form)
   Ikl = I(form,:);

   x = 1:winL;

   for Ndots = winL:-1:1
      [srt,ix] = sort(Ikl,'descend');

      % Guessing
      win = sort(ix(1:Ndots));

      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      for per = 1:perN
         sig(per,1:Ndots) = (nrm(f(per,win) - mean(f(per,win))));
         for btype = 1:btypeN
            cor(btype,per) = sig(per,1:Ndots) * port{btype}(win)';
            cor(btype,per) = (cor(btype,per) +1)/2;
         end

         [~,ind] = max(cor(:,per));
         if ~any(per == Bnum{3})
            des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
         end
      end
      des_m(form,winL - Ndots +1) = (des(1,1)+des(2,2))/2;
   end
end
%% Kullback
figure('Color','w')
plot(des_m(1,:),'.-k'),grid,axis tight,ylim([.5 1])
xlabel('Quantity of excluded samples','FontName','Times New Roman','FontSize',14)
ylabel('Recognition quality','FontName','Times New Roman','FontSize',14)
%% Az
figure('Color','w')
plot(des_m(2,:),'.-k'),grid,axis tight,ylim([.5 1])
% xlabel('Quantity of excluded samples','FontName','Times New Roman','FontSize',14)
% ylabel('Recognition quality','FontName','Times New Roman','FontSize',14)
xlabel('Количество исключённых отсчётов','FontName','Times New Roman','FontSize',14)
ylabel('Качество распознавания','FontName','Times New Roman','FontSize',14)












