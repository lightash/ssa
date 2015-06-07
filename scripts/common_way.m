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
for per = 1:perN

   period = mark(per);
   window = period+Bwin(1): period+Bwin(2);
   f(per,:) = in(window);

   f(per,:) = nrm(f(per,:),1);
end

port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(Bnum{btype},:),1));
end

% Guessing
disp('Guessing')
% btypeN = 2;

des = zeros(btypeN);
for per = 1:perN
   disp(per)
   for btype = 1:btypeN
      cor(btype,per) = f(per,:) * port{btype}';
      cor(btype,per) = (cor(btype,per) +1)/2;
   end
   [~,ind] = max(cor(:,per));
      
%       [~,q12(per,:)] = des_MOD(port{1},port{2},f(per,:));
%       [~,q23(per,:)] = des_MOD(port{2},port{3},f(per,:));
%       [~,q13(per,:)] = des_MOD(port{1},port{3},f(per,:));
%       [~,ind] = max(q12(per,:));
      
%    if ~any(per == Bnum{3})
      des(Bord(per),ind) = des(Bord(per),ind) + 1/Blen(Bord(per));
%    end
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






