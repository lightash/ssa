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
port{4} = nrm(mean([port{1};port{2};port{3}],1));
%%
% Detecting
disp('Detecting')
for btype = 1:btypeN+1
   disp(btype)
   for pos = 1:length(in)-winL+1
      win = pos:pos+winL-1;
      
      sig = in(win);
      sig = nrm(sig - mean(sig));
      
      cor(btype,pos) = port{btype} * sig';
      cor(btype,pos) = (cor(btype,pos) +1)/2;
   end
end
%%
plot(cor'),hold on
plot(mark-46,.8*ones(1,length(mark)),'.')





















