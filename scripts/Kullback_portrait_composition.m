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
win = [-47 80];   % Borders of PQRST period
winL = win(2)-win(1)+1;

% Generating portraits
disp('Generating portraits')
perN = all_beats(end);  % Number of periods to use
f = zeros(perN, winL );
for per = 1:perN

   period = mark(per);
   window = period+win(1): period+win(2);
   f(per,:) = in(window);
   f(per,:) = f(per,:) - mean(f(per,:));

   f(per,:) = nrm(f(per,:));
end
port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(Bnum{btype},:),1));
end

btypeN = 2;
% %%
% for nb = 1:8
Nbins = 2;%^nb;
% clear hyp scale H
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
% Ikl = DK( H{1}, H{2}, Blen(1), Blen(2) );
Ikl = AlphaZ( H{1}, H{2} );

load('indei_NA.mat')
[~,maxi] = max(indei);
% x = 1:winL;
% figure,plot(x,port{1},'-',x,port{2},'-g')%,x,port{3},'-r')
% hold on,plot(wini{maxi-1},port{1}(wini{maxi-1}),'.',wini{maxi-1},port{2}(wini{maxi-1}),'.g')%wini{maxi-1},port{3}(wini{maxi-1}),'.r')
% grid,axis tight
% figure,subplot(211),plot(indei,'.-'),axis tight,subplot(212),plot(maxj,'.-'),axis tight

Ndots = length(wini{maxi-1});
[srt,ix] = sort(Ikl,'descend');
% hold on,plot(x(sort(ix(1:Ndots))),port{1}(sort(ix(1:Ndots))),'o')

% Guessing
disp('Guessing')
win = sort(ix(1:Ndots))';

cor = zeros(btypeN,perN);
des = zeros(btypeN);
for per = 1:perN
   disp(per)
   sig(per,:) = (nrm(f(per,win) - mean(f(per,win))));
   for btype = 1:btypeN
      cor(btype,per) = sig(per,:) * port{btype}(win)';
      cor(btype,per) = (cor(btype,per) +1)/2;
   end
   
   [~,ind] = max(cor(:,per));
   if ~any(per == Bnum{3})
      des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
   end
end
% binn(nb) = (des(1,1)+des(2,2))/2;
% end
% %%
figure
k = 0;
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
      subplot(btypeN,btypeN,k),stem(des(i,j),'.-'),axis([0 2 0 1])
      xlabel(des(i,j))
  end
end
title((des(1,1)+des(2,2))/2)



