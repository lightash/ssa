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
   f(per,:) = f(per,:) - mean(f(per,:));

   f(per,:) = nrm(f(per,:));
end
port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(Bnum{btype},:),1));
end

%% pic 2 Orthogonal components
dec = MOD(port{1},port{2});

x1 = 1:winL;
x2 = winL+1:2*winL;
x3 = 2*winL+1:3*winL;
figure('color','white')
plot(x1,port{1},'k','LineWidth',1.5),hold on
plot(x1,dec.s21*1.5,'--k','LineWidth',1.5)
plot(x2,port{2},'k','LineWidth',1.5)
plot(x2,dec.s12*1.5,'--k','LineWidth',1.5)
plot(x3,dec.s11*1.5,'k','LineWidth',1.5)
plot(x3,dec.s22*1.5,'--k','LineWidth',1.5)
plot([x1(end) x1(end)],[-.3 .7],'--k','LineWidth',3)
plot([x2(end) x2(end)],[-.3 .7],'--k','LineWidth',3),grid,axis tight
leg = legend('S_1','s_2_1','S_2','s_1_2','s_1_1','s_2_2');
set(leg,'FontName','Times New Roman','FontSize',14)
xlabel({'Сквозная нумерация отсчётов';'Коллинеарные составляющие 1     Коллинеарные составляющие 2     Эксклюзивные составляющие'},...
   'FontName','Times New Roman','FontSize',14)
ylabel('Ортогональные составляющие сигнала','FontName','Times New Roman','FontSize',14)

%% pic 3 Histograms samples
for per = 1:perN
   dec = MOD(port{1},port{2},f(per,:));
   OrthDec{1}(per,:) = dec.x1;
   OrthDec{2}(per,:) = dec.x2;
   OrthDec{3}(per,:) = dec.x11;
   OrthDec{4}(per,:) = dec.x22;
end

btypeN = 2;
Nbins = 16;
% for collinear x1 and x2
for j = 1:2
   for i = 1:winL
      [hyp(i,:),scale(i,:)] = hist( [OrthDec{j}(Bnum{1},i);OrthDec{j}(Bnum{2},i)] ,Nbins);

      for btype = 1:btypeN
         H{j,btype}(:,i) = hist( OrthDec{j}(Bnum{btype},i), scale(i,:) )/Blen(btype)';
         % Распределить наименьший столбец по нулям
         [minh,indh] = min(H{j,btype}( H{j,btype}(:,i)>0, i));
         mzh = minh/( length( H{j,btype}( H{j,btype}(:,i)==0, i)) +1);
         H{j,btype}( H{j,btype}(:,i)==0 | H{j,btype}(:,i)==minh, i) = mzh;
      end
   end
end
% for exclusive x11 and x22
for j = 3
   for i = 1:winL
      [hyp(i,:),scale(i,:)] = hist( [OrthDec{3}(Bnum{1},i);OrthDec{4}(Bnum{2},i)] ,Nbins);

      for btype = 1:btypeN
         H{j,btype}(:,i) = hist( OrthDec{btype+2}(Bnum{btype},i), scale(i,:) )/Blen(btype)';
         % Распределить наименьший столбец по нулям
         [minh,indh] = min(H{j,btype}( H{j,btype}(:,i)>0, i));
         mzh = minh/( length( H{j,btype}( H{j,btype}(:,i)==0, i)) +1);
         H{j,btype}( H{j,btype}(:,i)==0 | H{j,btype}(:,i)==minh, i) = mzh;
      end
   end
end

%% pic 3 graph
figure('Color','w')

i = 20;  % Hist to display
for k = 1:3
   subplot(1,3,k)
   
   Hb = zeros(4,size(H{k,1}(:,i),1)); % Nback Aback Nfore Afore
   for j = 1:size(H{k,1}(:,i),1)
      if H{k,1}(j,i) > H{k,2}(j,i)
         Hb(1,j) = H{k,1}(j,i);
         Hb(4,j) = H{k,2}(j,i);
      else
         Hb(3,j) = H{k,1}(j,i);
         Hb(2,j) = H{k,2}(j,i);
      end
   end
   
   obj = bar(scale(i,:),Hb(1,:));hold on
   set(obj,'FaceColor',[.5 .5 .5]);
   bar(scale(i,:),Hb(2,:),'w')
   obj = bar(scale(i,:),Hb(3,:));
   set(obj,'FaceColor',[.5 .5 .5]);
   bar(scale(i,:),Hb(4,:),'w')
   obj = legend('N','A','Location','North');
   set(obj,'FontName','Times New Roman','FontSize',14)
   xlabel('Интервалы значений составляющих','FontName','Times New Roman','FontSize',14)
   if k == 1
      ylabel({'Relative frequency for' 'the intervals of values'},'FontName','Times New Roman','FontSize',14)
   end
   grid,axis tight,ylim([0 .3])
end

%% pic 4 Kullback informativity
for i = 1:size(H,1)
   I{1}(i,:) = DK( H{i,1}, H{i,2}, Blen(1), Blen(2) );
end

figure('Color','w')
plot([I{1}(1,1) I{1}(2,1) I{1}(3,:)],'.-k'),grid,axis tight
title('Kullback')
xlabel('Номер отсчёта','FontName','Times New Roman','FontSize',14)
ylabel('Информативность','FontName','Times New Roman','FontSize',14)

%% pic 5 Alpha z informativity
for i = 1:size(H,1)
   I{2}(i,:) = AlphaZ( H{i,1}, H{i,2} );
end

figure('Color','w')
subplot(211)
plot([I{2}(1,1) I{2}(2,1) I{2}(3,:)],'.-k'),grid,axis tight,ylim([.46 .49])
title('\alpha_z')
xlabel('Номер отсчёта','FontName','Times New Roman','FontSize',14)
ylabel('Информативность','FontName','Times New Roman','FontSize',14)
subplot(212)
plot([I{2}(1,1) I{2}(2,1) I{2}(3,:)],'.-k'),grid,axis tight
title('\alpha_z')
xlabel('Номер отсчёта','FontName','Times New Roman','FontSize',14)
ylabel('Информативность','FontName','Times New Roman','FontSize',14)

%% pic 6 Exclusion informativity, first iteration, excluding after orthogonal decomposition
NAbeats = [Bnum{1} Bnum{2}];
NAbeatsN = Blen(1)+Blen(2);
dec = cell(1,NAbeatsN);
for per = 1:NAbeatsN
   dec{per} = MOD(port{1},port{2},f(NAbeats(per),:));
end

% Full probability
des = zeros(btypeN);
for per = 1:NAbeatsN
   ind = des_MOD(dec{per});
   des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
end
prob_first(:,1) = [des(1);des(end)];

init_win{1} = 1:winL;
init_win{2} = 1:winL;
init_win{3} = 1:winL;

% Excuding by one
for exclude = 1:3*winL
   disp(exclude)
   
   win = init_win;
   
   part = ceil(exclude/winL);  % exclude belongs to collinear 1, 2 or exclusive
   local_excl = exclude - (part-1)*winL;
   
   win{part} = [1:local_excl-1 local_excl+1:winL];

   des = zeros(btypeN);
   for per = 1:NAbeatsN
      ind = des_MOD(dec{per},win);
      des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
   end
   prob_first(:,exclude+1) = [des(1);des(end)];
end
save('orth_inf_pic6','prob_first')

%% pic 6 graph
figure('color','white')
plot(prob_first(1,:),'k','LineWidth',1.5),hold on
plot(prob_first(2,:),'--k','LineWidth',1.5)
plot(mean(prob_first,1),':k','LineWidth',1.5)
plot([x1(end) x1(end)],[.5 1],'--k','LineWidth',3)
plot([x2(end) x2(end)],[.5 1],'--k','LineWidth',3),grid,axis tight
leg = legend('N','A','mean');
set(leg,'FontName','Times New Roman','FontSize',14)
xlabel({'Сквозная нумерация отсчётов';'Коллинеарные составляющие 1     Коллинеарные составляющие 2     Эксклюзивные составляющие'},...
   'FontName','Times New Roman','FontSize',14)
ylabel('Orth first iter probs','FontName','Times New Roman','FontSize',14)

%% pic 7 Delta of pic 6
dprob_first(1,:) = prob_first(1,:) - prob_first(1);
dprob_first(2,:) = prob_first(2,:) - prob_first(2);

figure('color','white')
plot(dprob_first(1,:),'k','LineWidth',1.5),hold on
plot(dprob_first(2,:),'-.k','LineWidth',1.5)
plot(mean(dprob_first,1),':k','LineWidth',1.5)
plot([x1(end) x1(end)],[-.5 .2],'--k','LineWidth',3)
plot([x2(end) x2(end)],[-.5 .2],'--k','LineWidth',3),grid,axis tight
leg = legend('N','A','mean');
set(leg,'FontName','Times New Roman','FontSize',14)
xlabel({'Сквозная нумерация отсчётов';'Коллинеарные составляющие 1     Коллинеарные составляющие 2     Эксклюзивные составляющие'},...
   'FontName','Times New Roman','FontSize',14)
ylabel('Orth first iter delta probs','FontName','Times New Roman','FontSize',14)

%% pic 8 sorted pic 7
figure('color','white')
plot(sort(dprob_first(1,:),'descend'),'k','LineWidth',1.5),hold on
plot(sort(dprob_first(2,:),'descend'),'-.k','LineWidth',1.5)
plot(sort(mean(dprob_first,1),'descend'),':k','LineWidth',1.5)
plot([x1(end) x1(end)],[-.5 .2],'--k','LineWidth',3)
plot([x2(end) x2(end)],[-.5 .2],'--k','LineWidth',3),grid,axis tight
leg = legend('N','A','mean');
set(leg,'FontName','Times New Roman','FontSize',14)
xlabel('Сквозная нумерация отсчётов','FontName','Times New Roman','FontSize',14)
ylabel('Orth first iter delta probs','FontName','Times New Roman','FontSize',14)

%% pic 9 Recognition quality by synchronous portraits composition, all iterations of exclusion informativity, excluding after orthogonal decomposition
% Full probability
des = zeros(btypeN);
for per = 1:NAbeatsN
   ind = des_MOD(dec{per});
   des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
end
prob_all(1) = mean([des(1) des(end)]);

% Excuding serially
win_rem = cell(1,winL-1);
win_rem{1} = 1:winL;
order = zeros(1,winL-1);
for throw = 1:winL-1
   
   prob_one = zeros(1,winL-throw+1);
   for exclude = 1:winL-throw+1
      disp([throw exclude])

      win = [win_rem{throw}(1:exclude-1) win_rem{throw}(exclude+1:end)];

      des = zeros(btypeN);
      for per = 1:NAbeatsN
         ind = des_MOD(MOD(nrm(port{1}(win)), nrm(port{2}(win)), f(NAbeats(per),win)));
         des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
      end
      prob_one(exclude) = mean([des(1) des(end)]);
   end
   
   [prob_all(throw),ind] = max(prob_one);
   win_rem{throw+1} = [win_rem{throw}(1:ind-1) win_rem{throw}(ind+1:end)];
   order(throw) = find( port{1} == port{1}(win_rem{throw}(ind)) );
end
order(winL) = win_rem{end};
save('orth_inf_pic9','prob_all','win_rem','order')

%% pic 9 graph
figure('color','white')
plot(prob_all,'.-k','LineWidth',1.5),grid,axis tight,ylim([.5 1])
xlabel('Сквозная нумерация отсчётов','FontName','Times New Roman','FontSize',14)
ylabel('Orth all iter mean probs','FontName','Times New Roman','FontSize',14)

%% pic 10 Recognition quality by asynchronous portraits composition, all iterations by first exclusion informativity (pic 6), excluding before orthogonal decomposition
load('orth_inf_pic6','prob_first')
[~,order] = sort(mean(prob_first(:,2:end),1),'descend');

% Full probability
des = zeros(btypeN);
for per = 1:NAbeatsN
   ind = des_MOD(dec{per});
   des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
end
prob_all(1) = mean([des(1) des(end)]);

% Excuding serially
win_rem = cell(winL-1,3);
win_rem{1,1} = 1:winL;
win_rem{1,2} = 1:winL;
win_rem{1,3} = 1:winL;
for throw = 1:3*winL-1
   disp(throw)
   
   part = ceil(order(throw)/winL);  % throw belongs to collinear 1, 2 or exclusive
   if part == 1
      bef = sum(order(1:throw-1) < order(throw));
   elseif part == 2
      bef = sum(order(1:throw-1) < order(throw) & order(1:throw-1) > winL);
   else
      bef = sum(order(1:throw-1) < order(throw) & order(1:throw-1) > 2*winL);
   end
   local_ind = order(throw) - (part-1)*winL - bef;

   win = win_rem(throw,:);
   win{part} = [win_rem{throw,part}(1:local_ind-1) win_rem{throw,part}(local_ind+1:end)];

   des = zeros(btypeN);
   for per = 1:NAbeatsN
      ind = des_MOD(dec{per},win);
      des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
   end
   
   prob_all(throw+1) = mean([des(1) des(end)]);
   
   win_rem(throw+1,:) = win;
end
save('orth_inf_pic10','prob_all','win_rem','order')

%% pic 10 graph
figure('color','white')
plot(prob_all,'.-k','LineWidth',1.5),grid,axis tight,ylim([.5 1])
xlabel('Сквозная нумерация отсчётов','FontName','Times New Roman','FontSize',14)
ylabel('Orth all serial iter mean probs','FontName','Times New Roman','FontSize',14)

%% pic - Recognition quality by asynchronous portraits composition, all iterations of exclusion informativity, excluding before orthogonal decomposition
% Full probability
des = zeros(btypeN);
for per = 1:NAbeatsN
   ind = des_MOD(dec{per});
   des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
end
prob_all(:,1) = mean([des(1) des(end)]);

% Excuding serially
win_rem = cell(winL-1,3);
win_rem{1,1} = 1:winL;
win_rem{1,2} = 1:winL;
win_rem{1,3} = 1:winL;
order = zeros(1,3*winL-1);
for throw = 1:3*winL-1
   disp(throw)
   tic
   prob_one = zeros(1,3*winL-throw+1);
   for exclude = 1:3*winL-throw+1
      
      win = win_rem(throw,:);

      if exclude <= size(win_rem{throw,1},2)  % collinear 1
         part = 1;
         local_excl = exclude;
      elseif exclude <= size(win_rem{throw,1},2) + size(win_rem{throw,2},2)  % collinear 2
         part = 2;
         local_excl = exclude - size(win_rem{throw,1},2);
      elseif exclude <= size(win_rem{throw,1},2) + size(win_rem{throw,2},2) + size(win_rem{throw,3},2)  % exclusive
         part = 3;
         local_excl = exclude - size(win_rem{throw,1},2) - size(win_rem{throw,2},2);
      end
      
      win{part} = [win_rem{throw,part}(1:local_excl-1) win_rem{throw,part}(local_excl+1:end)];
      
      des = zeros(btypeN);
      for per = 1:NAbeatsN
         ind = des_MOD(dec{per},win);
         des(Bord(NAbeats(per)),ind) = des(Bord(NAbeats(per)),ind) + 1/Blen(Bord(NAbeats(per)));
      end
      prob_one(local_excl) = mean([des(1) des(end)]);
   end
   
   [prob_all(throw),ind] = max(prob_one);
   
   if ind <= size(win_rem{throw,1},2)  % collinear 1
      part = 1;
      ind = ind;
   elseif ind <= size(win_rem{throw,1},2) + size(win_rem{throw,2},2)  % collinear 2
      part = 2;
      ind = ind - size(win_rem{throw,1},2);
   elseif ind <= size(win_rem{throw,1},2) + size(win_rem{throw,2},2) + size(win_rem{throw,3},2)  % exclusive
      part = 3;
      ind = ind - size(win_rem{throw,1},2) - size(win_rem{throw,2},2);
   end
   
   win_rem(throw+1,:) = win_rem(throw,:);
   win_rem{throw+1,part} = [win_rem{throw,part}(1:ind-1) win_rem{throw,part}(ind+1:end)];
   
   order(throw) = win_rem{throw,part}(ind) + (part-1)*winL;
end
order(3*winL) = (part-1)*winL + win_rem{end,part}(ind);
save('orth_inf_pic-','prob_all','win_rem','order')
% system('shutdown -h')

%% pic - graph
figure('color','white')
plot(prob_all,'.-k','LineWidth',1.5),grid,axis tight,ylim([.5 1])
xlabel('Сквозная нумерация отсчётов','FontName','Times New Roman','FontSize',14)
title('Recognition quality by asynchronous portraits composition, all iterations of exclusion informativity, excluding before orthogonal decomposition','FontName','Times New Roman','FontSize',14)













