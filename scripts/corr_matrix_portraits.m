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

btypeN = 2;
cormatL = (winL^2-winL)/2;

% Generating portraits
disp('Generating portraits')
perN = all_beats(end);  % Number of periods to use
f = zeros(perN, winL );
fcm = zeros(perN,cormatL);
for per = 1:perN
   disp(per)
   period = mark(per);
   window = period+win(1): period+win(2);
   f(per,:) = in(window);
   f(per,:) = f(per,:) - mean(f(per,:));
   f(per,:) = nrm(f(per,:));
   
   k = 0;
   for i = 1:winL
      for j = i+1:winL
         k = k+1;
         fcm(per,k) = f(per,i) * f(per,j);
      end
   end
   fcm(per,:) = nrm(fcm(per,:) - mean(fcm(per,:)));
end
port = cell(1,btypeN);
for btype = 1:btypeN
   nmf = nrm(mean(f(Bnum{btype},:),1));
   port{btype} = zeros(1,(winL^2-winL)/2);
   k = 0;
   for i = 1:winL
      for j = i+1:winL
         k = k+1;
         port{btype}(k) = nmf(i) * nmf(j);
      end
   end
   port{btype} = nrm(port{btype} - mean(port{btype}));
   
%    port{btype} = mean(fcm(Bnum{btype},:),1);
%    port{btype} = nrm(port{btype} - mean(port{btype}));
end
%%
% Composition
for nb = 1:6
   disp(nb)
Nbins = 2^nb;%16;
clear hyp scale H
for i = 1:cormatL
   [hyp(i,:),scale(i,:)] = hist( [fcm(Bnum{1},i);fcm(Bnum{2},i)] ,Nbins);
   
   for btype = 1:btypeN
      H{btype}(:,i) = hist( fcm(Bnum{btype},i) ,scale(i,:))/Blen(btype)';
      % Распределить наименьший столбец по нулям
      [minh,indh] = min(H{btype}( H{btype}(:,i)>0, i));
      mzh = minh/( length( H{btype}( H{btype}(:,i)==0, i)) +1);
      H{btype}( H{btype}(:,i)==0 | H{btype}(:,i)==minh, i) = mzh;
   end
end
inf_dk = DK( H{1}, H{2}, Blen(1), Blen(2) );
inf_az = AlphaZ( H{1}, H{2} );
% %%
[srt(nb,:),isrt(nb,:)] = sort(inf_dk,'descend');
[srt1(nb,:),isrt1(nb,:)] = sort(inf_az,'descend');
end
%%
nd_max = 1500;
sN = zeros(1,nd_max);
for nd = 1:nd_max
   Ndots = 1*nd;
   ixN = isrt(:,1:Ndots);

   ixL = zeros(2,cormatL);
   ixL(1,ixN(1,:)) = 1;
   ixL(2,ixN(2,:)) = 1;

   ixA = ixL(1,:) & ixL(2,:);
   sN(nd) = sum(ixA)/Ndots;
end
% plot(sN)
% figure
% plot(port{1}),hold on
% plot(port{2},'g'),hold on
% plot(ixN(1,:),port{1}(ixN(1,:)),'.'),hold on
% plot(ixN(2,:),port{1}(ixN(2,:)),'og'),axis tight
% hold on
% plot(ixL'),hold on
% plot(ixA,'r')

%%
% Guessing
disp('Guessing')

des = zeros(btypeN,btypeN,nd_max);
for nd = 1:nd_max
   disp(nd)
   Ndots = 1*nd;
   win = sort(isrt(2,1:Ndots));
   for per = 1:perN
   %    disp(per)

         for btype = 1:btypeN
            cor(btype,per) = fcm(per,win) * port{btype}(win)';
            cor(btype,per) = (cor(btype,per) +1)/2;
         end

      if ~any(per == Bnum{3})
         [~,ind] = max(cor(:,per));
         des(Bord(per),ind,nd) = des(Bord(per),ind,nd) + 1/Blen(Bord(per));
      end
   end
end
%%
figure
k = 0;
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
      subplot(btypeN,btypeN,k),plot( permute(des(i,j,:),[1 3 2]) ,'.-'),axis tight%([0 2 0 1])
%       xlabel(mean(temp))
  end
end
% title((des(1,1)+des(2,2)+des(3,3))/3)
% title((des(1,1)+des(2,2))/2)
figure,plot(permute((des(1,1,:)+des(2,2,:))/2,[1 3 2]))
% %%
figure
plot(port{1}),hold on
plot(win,port{1}(win),'.'),hold on
plot(port{2},'g'),hold on
plot(win,port{2}(win),'.g'),hold on
% plot(port{3},'r')
axis tight
% figure
% subplot(321),plot(port{1}),axis tight
% subplot(323),plot(port{2},'g'),axis tight
% subplot(325),plot(port{3},'r'),axis tight
% subplot(322),plot(cms(Bnum{1},:)'),axis tight
% subplot(324),plot(cms(Bnum{2},:)'),axis tight
% subplot(326),plot(cms(Bnum{3},:)'),axis tight











