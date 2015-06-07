clc
% close all
clear all

load('D:\Dropbox\Signals\incartdb\I20\I20proc.mat')
in = val(1,:);
annot(2461) = 'N';   % Fusion of ventricular and normal beat in 'I20'

all_beats = 1:length(annot);
% Normal beat, Atrial premature beat, Premature ventricular contraction
bmark = 'NAV';
btypeN = 3;       % Beat types to examine
for i = 1:btypeN
   Bnum{i} = all_beats(annot == bmark(i));
   Blen(i) = length(Bnum{i});
   Bpos{i} = mark(Bnum{i});
   Bord(all_beats(annot ==  bmark(i))) = i*ones(1,length(Bpos{i}));
end
Bwin = [-47 80];   % Borders of PQRST
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
% port = cell(1,btypeN);
% for btype = 1:btypeN
% %    port{btype} = nrm(mean(f(Bnum{btype},:),1));
%    port{btype} = nrm(AM(f(Bnum{btype},:)),1);
% end

%%
% [e, c, a2] = AM(f);
% ec = (e'*c)';
% ecv = zeros(size(in));
% fv = zeros(size(in));
% pv = zeros(size(in));
% a2v = zeros(size(in));
% for per = 1:perN
%    ecv(mark(per)+Bwin(1):mark(per)+Bwin(2)) = nrm(ec(per,:));
%    fv(mark(per)+Bwin(1):mark(per)+Bwin(2)) = f(per,:);
%    pv(mark(per)+Bwin(1):mark(per)+Bwin(2)) = port{1};
%    a2v(mark(per)+Bwin(1):mark(per)+Bwin(2)) = a2(per,:);
% end
% figure
% plot(fv),hold on
% plot(pv,'r'),hold on
% plot(ecv,'g'),hold on
% plot(a2v,':k'),hold on
% plot(ecv+a2v,'y'),hold on
%%
[E, C, A2] = impAM(f(Bnum{1}(1:winL),:),'from_end');
%%
figure,plot(E')















