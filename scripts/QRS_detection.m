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
   f(per,:) = f(per,:) - mean(f(per,:));

   f(per,:) = nrm(f(per,:));
end

port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(Bnum{btype},:),1));
end
%%
% Detecting
disp('Detecting')
cor = zeros(btypeN,length(in));
for btype = 1:btypeN
   disp(btype)
   for pos = 1:length(in)-winL+1
      win = pos:pos+winL-1;
      
      sig = in(win);
      sig = nrm(sig - mean(sig));
      
      cor(btype,pos-Bwin(1)) = port{btype} * sig';
      cor(btype,pos-Bwin(1)) = (cor(btype,pos-Bwin(1)) +1)/2;
   end
end
%%
figure
plot(cor'),hold on
plot(mark,.8*ones(1,length(mark)),'.'),hold on
plot(in+2,'m')
%%
% % Global threshold N27 A192 V162
% inx = 1:length(in);
% thrN = 2;
% gthr = linspace(.7,.97,thrN);
% for thri = 1:thrN
%    for btype = 1:btypeN
%       det(thri,btype) = sum(cor(btype,:)>=gthr(thri));
%       cord{thri,btype} = inx(cor(btype,:)>=gthr(thri));
%    end
% end
%%
% Local maximum
lm_win = round(1.3*mean(diff(mark)));

for btype = 1%:btypeN
   disp(btype)

   cord = zeros(1,length(cor(btype,:))-lm_win+1);
   for pos = 1:length(cor(btype,:))-lm_win+1
      win = pos:pos+lm_win-1;
      
      [~,mi] = max(cor(btype,win));
      cord(pos) = pos+mi-1;
   end
   
   cord = unique(cord);
end

%%
x = [gthr;gthr;gthr]';
figure,plot(x,det)
title(Blen'),axis tight
%%
figure
stem(cord,1.02*ones(1,length(cord)),'m'),hold on
stem(Bpos{1},1.01*ones(1,length(Bpos{1})),'.k'),hold on
stem(Bpos{2},1.01*ones(1,length(Bpos{2})),'sk'),hold on
stem(Bpos{3},1.01*ones(1,length(Bpos{3})),'dk'),hold on
% stem(cord{27,1},ones(1,length(cord{27,1})),'b'),hold on
% stem(cord{192,1},.99*ones(1,length(cord{192,1})),'g'),hold on
% stem(cord{162,1},.98*ones(1,length(cord{162,1})),'r'),hold on
ylim([.9 1.1])











