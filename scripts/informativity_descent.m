clc;
close all;
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

%% Informativity descent
cor = zeros(btypeN,perN);
des = zeros(btypeN);
indei = zeros(1,winL);
maxj = zeros(1,winL);
wini = cell(1,winL);
wini{1} = 1:winL;
fg = f;

for per = 1:perN
   for btype = 1:btypeN
      cor(btype,per) = fg(per,:) * port{btype}';
   end
   [~,ind] = max(cor(:,per));
   des(Bord(per),ind) = des(Bord(per),ind) + 1/Blen(Bord(per));
end
indei(1) = (des(1,1)+des(2,2)+des(3,3))/3;

Nbins = 100;
for i = 1:winL
%    H{btype} = hist(fg(Bnum{btype},:),Nbins)/Blen(btype);
   T = hist( [fg(Bnum{1},i) fg(Bnum{2},i) fg(Bnum{3},i)] ,Nbins);
end
Ikl(:,1:2) = DKL( H{1}, H{2} );
Ikl(:,3:4) = DKL( H{1}, H{3} );
Ikl(:,5:6) = DKL( H{2}, H{3} );
%%
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
         des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
         
      end
      indej(j) = (des(1,1)+des(2,2)+des(3,3))/3;
%       DKL{i} = DKL(fg(:,winj));
      
   end
   
   [indei(i),maxj(i)] = max(indej);

   wini{i} = [wini{i-1}(1:maxj(i)-1) wini{i-1}(maxj(i)+1:end)];
      
end

x = 1:winL;
[~,maxi] = max(indei);
figure,plot(x,port{1},'-',x,port{2},'-g',x,port{3},'-r')
hold on,plot(wini{maxi-1},port{1}(wini{maxi-1}),'.',wini{maxi-1},port{2}(wini{maxi-1}),'.g',...
   wini{maxi-1},port{3}(wini{maxi-1}),'.r'),grid,axis tight
figure,subplot(211),plot(indei,'.-'),axis tight,subplot(212),plot(maxj,'.-'),axis tight