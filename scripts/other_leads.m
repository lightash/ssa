clc;
clear all;

file = '21';
load(['D:\Dropbox\Signals\incartdb\I' file '\I' file 'proc.mat'])
% annot(2461) = 'N';   % Fusion of ventricular and normal beat

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
perN = all_beats(end);  % Number of periods to use

figure(1)
des = zeros(size(val,1),btypeN,btypeN);
for lead = 1:size(val,1)
   in = val(lead,:);

   % Generating portraits

   f = zeros(perN, winL );
   for per = 1:perN
      window = mark(per)+Bwin(1): mark(per)+Bwin(2);
      f(per,:) = in(window);
      f(per,:) = nrm(f(per,:),1);
   end

   port = cell(1,btypeN);
   for btype = 1:btypeN
      port{btype} = nrm(mean(f(Bnum{btype},:),1));
   end

   x = (lead-1)*3*winL+1:lead*3*winL;
   plot(x,[port{1} port{2} port{3}]+.5),hold on

   % Guessing
   for per = 1:perN

      for btype = 1:btypeN
         cor(btype,per) = (f(per,:) * port{btype}' +1)/2;
      end

      [~,ind] = max(cor(:,per));
      des(lead,Bord(per),ind) = des(lead,Bord(per),ind) + 1/Blen(Bord(per));
   end

end
grid,axis tight
%%
figure
k = 0;
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
      subplot(btypeN,btypeN,k),stem(des(:,i,j),'.-'),axis tight,ylim([0 1])
      xlabel(des(i,j))
   end
end
title((des(1,1)+des(2,2)+des(3,3))/3)











