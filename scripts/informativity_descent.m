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
probi = zeros(btypeN,2*winL);
maxj = zeros(1,2*winL);
order = zeros(2,2*winL);
wini = cell(btypeN,2*winL);
wini{1,1} = 1:winL;
wini{2,1} = 1:winL;

for per = 1:perN
   for btype = 1:btypeN
      cor(btype,per) = f(per,:) * port{btype}';
   end
   [~,ind] = max(cor(:,per));
   if ~any(per == Bnum{3})
      des(Bord(per),ind) = des(Bord(per),ind) + 1/Blen(Bord(per));
%       des(Bord(per),ind) = des(Bord(per),ind) + 1/(Blen(1)+Blen(2));
   end
end
% probi(1) = (des(1,1)+des(2,2)+des(3,3))/btypeN;
probi(1,1) = des(1,1);
probi(2,1) = des(2,2);
for i = 2:2*winL
   disp(i)
   tic
   
   probj{1} = zeros(1,size(wini{1,i-1},2));
   probj{2} = zeros(1,size(wini{2,i-1},2));
   for j = 1:size(wini{1,i-1},2)  % In N

      winj = [wini{1,i-1}(1:j-1) wini{1,i-1}(j+1:end)];

      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      for per = 1:perN

         cor(1,per) = nrm(f(per,winj),1) * nrm(port{btype}(winj),1)';
         cor(2,per) = nrm(f(per,wini{2,i-1}),1) * nrm(port{btype}(wini{2,i-1}),1)';

         [~,ind] = max(cor(:,per));
         if ~any(per == Bnum{3})
            des(Bord(per),ind) = des(Bord(per),ind) + 1/Blen(Bord(per));
%             des(Bord(per),ind) = des(Bord(per),ind) + 1/(Blen(1)+Blen(2));
         end
      end
      probj{1}(j) = des(1,1);
      probj{2}(j) = des(2,2);
   end
   
   for j = 1:size(wini{2,i-1},2)  % In A

      winj = [wini{2,i-1}(1:j-1) wini{2,i-1}(j+1:end)];

      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      for per = 1:perN

         cor(1,per) = nrm(f(per,wini{1,i-1}),1) * nrm(port{btype}(wini{1,i-1}),1)';
         cor(2,per) = nrm(f(per,winj),1) * nrm(port{btype}(winj),1)';

         [~,ind] = max(cor(:,per));
         if ~any(per == Bnum{3})
            des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
%             des(Bord(per),ind) = des(Bord(per),ind) + 1/(Blen(1)+Blen(2));
         end
      end
      probj{1}(size(wini{1,i-1},2) + j) = des(1,1);
      probj{2}(size(wini{1,i-1},2) + j) = des(2,2);
   end
   
   dprobj_mean = mean( [probj{1}-probi(1,i-1) ; probj{2}-probi(2,i-1)] ,1);
   if size(wini{1,i-1},2) > 1
      if size(wini{2,i-1},2) > 1
         [~,maxj(i)] = max(dprobj_mean);
      else
         dprobj_mean(size(wini{1,i-1},2) + 1) = -1;
         [~,maxj(i)] = max(dprobj_mean);
      end
   else
      if size(wini{2,i-1},2) > 1
         dprobj_mean(1) = -1;
         [~,maxj(i)] = max(dprobj_mean);
      else
         [~,maxj(i)] = max(dprobj_mean);
      end
   end
   
   if maxj(i) <= size(wini{1,i-1},2)
      probi(1,i) = probj{1}(maxj(i));
      probi(2,i) = probj{2}(maxj(i));
      wini{1,i} = [wini{1,i-1}(1:maxj(i)-1) wini{1,i-1}(maxj(i)+1:end)];
      wini{2,i} = wini{2,i-1};
      order(1,i-1) = find(port{1} == port{1}(wini{1,i-1}(maxj(i))));
      order(2,i-1) = 1;
   else
      maxj(i) = maxj(i) - size(wini{1,i-1},2);
      probi(1,i) = probj{1}(maxj(i));
      probi(2,i) = probj{2}(maxj(i));
      wini{1,i} = wini{1,i-1};
      wini{2,i} = [wini{2,i-1}(1:maxj(i)-1) wini{2,i-1}(maxj(i)+1:end)];
      order(1,i-1) = find(port{1} == port{1}(wini{2,i-1}(maxj(i))));
      order(2,i-1) = 2;
   end
   
   t(i-1) = toc;
end
for i = 1:2*winL  % Filling the last one
   if size(find(order(1,:)==i),2)==0  % If no such 'i' in 'order'
      if i <= winL
         order(1,2*winL) = i;
         order(2,2*winL) = 1;
      else
         order(1,2*winL) = i-winL;
         order(2,2*winL) = 2;
      end
   end
end
%%
probi_mean = mean(probi,1);
inform = zeros(2,winL);
for i = 1:2*winL-1
   inform(order(2,i),order(1,i)) = probi_mean(i+1) - probi_mean(i);
end
inform(order(2,2*winL),order(1,2*winL)) = .5 - probi_mean(2*winL);
% %%
% load('indei.mat')
% x = 1:winL;
% [~,maxi] = max(indei);
% figure,plot(x,port{1},'-',x,port{2},'-g')%,x,port{3},'-r')
% hold on,plot(wini{maxi-1},port{1}(wini{maxi-1}),'.',wini{maxi-1},port{2}(wini{maxi-1}),'.g')%wini{maxi-1},port{3}(wini{maxi-1}),'.r')
% grid,axis tight
% figure,subplot(211),plot(indei,'.-'),axis tight,subplot(212),plot(maxj,'.-'),axis tight
% %%
% % Guessing
% disp('Guessing')
% win = wini{maxi-1};
% % win = 1:winL;
% des = zeros(btypeN);
% f = f;
% for per = 1:perN
%    disp(per)
%       sig(per,:) = (nrm(f(per,win) - mean(f(per,win))));
%       for btype = 1:btypeN
%          cor(btype,per) = sig(per,:) * port{btype}(win)';
%          cor(btype,per) = (cor(btype,per) +1)/2;
%       end
%       
%       [~,ind] = max(cor(:,per));
%       if ~any(per == Bnum{3})
%          des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
%       end
% end
% figure
% k = 0;
% for i = 1:btypeN
%    for j = 1:btypeN
%       k = k+1;
%       
%       subplot(btypeN,btypeN,k),stem(des(i,j),'.-'),axis([0 2 0 1])
%       xlabel(des(i,j))
%   end
% end
% title((des(1,1)+des(2,2))/2)







