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
%% Informativity descent
cor = zeros(btypeN,perN);
des = zeros(btypeN);
indei = zeros(1,winL);
maxj = zeros(1,winL);
order = zeros(1,winL);
wini = cell(1,winL);
wini{1} = 1:winL;
fg = f;

for per = 1:perN
   for btype = 1:btypeN
      cor(btype,per) = fg(per,:) * port{btype}';
   end
   [~,ind] = max(cor(:,per));
   if ~any(per == Bnum{3})
      des(Bord(per),ind) = des(Bord(per),ind) + 1/Blen(Bord(per));
   end
end
% indei(1) = (des(1,1)+des(2,2)+des(3,3))/btypeN;
indei(1) = (des(1,1)+des(2,2))/btypeN;

for i = 2:winL
   disp(i)

   indej = zeros(1,length(wini{i-1}));
   for j = 1:length(wini{i-1})

      winj = [wini{i-1}(1:j-1) wini{i-1}(j+1:end)];

      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      for per = 1:perN
         
         for btype = 1:btypeN
            cor(btype,per) = nrm(fg(per,winj) - mean(fg(per,winj))) * ...
               nrm(port{btype}(winj) - mean(port{btype}(winj)))';
         end
         
         [~,ind] = max(cor(:,per));
         if ~any(per == Bnum{3})
            des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
         end
      end
%       indej(j) = (des(1,1)+des(2,2)+des(3,3))/btypeN;
      indej(j) = (des(1,1)+des(2,2))/btypeN;
%       DKL{i} = DKL(fg(:,winj));
      
   end
   
   [indei(i),maxj(i)] = max(indej);
   wini{i} = [wini{i-1}(1:maxj(i)-1) wini{i-1}(maxj(i)+1:end)];
   order(i-1) = find(port{1} == port{1}(wini{i-1}(maxj(i))));
end
for i = 1:winL
   if size(find(order==i),2)==0
      order(winL)=i;
   end
end
%%
load('indei_NA.mat')
x = 1:winL;
[~,maxi] = max(indei);
% figure,plot(x,port{1},'-',x,port{2},'-g')%,x,port{3},'-r')
% hold on,plot(wini{maxi-1},port{1}(wini{maxi-1}),'.',wini{maxi-1},port{2}(wini{maxi-1}),'.g')%wini{maxi-1},port{3}(wini{maxi-1}),'.r')
% grid,axis tight
% figure,subplot(211),plot(indei,'.-'),axis tight,subplot(212),plot(maxj,'.-'),axis tight
%%
% Guessing
disp('Guessing')
win = wini{maxi-1};
% win = 1:winL;
des = zeros(btypeN);
fg = f;
for per = 1:perN
   disp(per)
      sig(per,:) = (nrm(fg(per,win) - mean(fg(per,win))));
      for btype = 1:btypeN
         cor(btype,per) = sig(per,:) * port{btype}(win)';
         cor(btype,per) = (cor(btype,per) +1)/2;
      end
      
      [~,ind] = max(cor(:,per));
      if ~any(per == Bnum{3})
         des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
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
title((des(1,1)+des(2,2))/2)







