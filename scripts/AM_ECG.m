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
   f(per,:) = nrm(f(per,:),1);
end
btypeN = 2;
port = cell(1,btypeN);
for btype = 1:btypeN
%    port{btype} = nrm(mean(f(Bnum{btype},:),1));
   port{btype} = nrm(AM(f(Bnum{btype},:)),1);
end

% %%
% [e, c, a2] = AM(f);
% ec = c'*e;
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
E = cell(1,btypeN);
C = cell(1,btypeN);
Bas = cell(1,btypeN);
Nport = cell(1,btypeN);
for btype = 1:btypeN
   [E{btype}, C{btype}] = impAM(f(Bnum{btype},:),'from_end');
   [~,Bas{btype}] = GSOrth(E{btype});
   figure,plot(Bas{btype}')
   Nport{btype} = nrm((Bas{btype} * port{btype}')');
   figure,plot(Nport{btype})
end
%%
% Guessing
des = zeros(btypeN);
for per = 1:perN
   for btype = 1:btypeN
      Nf{btype}(per,:) = nrm((Bas{btype} * f(per,:)')');
      cor(btype,per) = Nf{btype}(per,:) * Nport{btype}';
      cor(btype,per) = (cor(btype,per) +1)/2;
   end
   [~,ind] = max(cor(:,per));
      
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
title((des(1,1)+des(2,2))/2)















