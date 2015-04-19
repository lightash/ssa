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
% figure
% plot(port{1}),hold on
% plot(port{2},'g')

btypeN = 2;

Nbins = 16;
for i = 1:winL
   [hyp(i,:),scale(i,:)] = hist( [f(Bnum{1},i);f(Bnum{2},i)] ,Nbins);
%    scalew(i) = (scale(i,2) - scale(i,1))/2;
   
   for btype = 1:btypeN
      H{btype}(:,i) = hist( f(Bnum{btype},i) ,scale(i,:))/Blen(btype)';
      % Распределить наименьший столбец по нулям
      [minh,indh] = min(H{btype}( H{btype}(:,i)>0, i));
      mzh = minh/( length( H{btype}( H{btype}(:,i)==0, i)) +1);
      H{btype}( H{btype}(:,i)==0 | H{btype}(:,i)==minh, i) = mzh;
   end
end

% Probabilities by histograms
for per = 1:perN
   if ~any(per == Bnum{3})
      for i = 1:winL
         [~,barnum] = min(abs(scale(i,:) - f(per,i)));
         hiprm(per,i) = H{Bord(per)}( barnum ,i);
      end
   end
end

for btype = 1:btypeN
%    hipr{btype} = mean(hiprm(Bnum{btype},:),2)';
   for i = 1:Blen(btype)
      [~,ix1] = sort(hiprm(Bnum{btype}(i),:),2,'descend');
      vol1 = sort(ix1(1:64));
      hipr{btype}(i) = mean(hiprm(Bnum{btype}(i),vol1),2)';
   end
end

% figure,plot(hipr)

figure
vom = 20;
for vo = 1:vom
   Ndots = round(vo*Blen/vom);
   for btype = 1:btypeN
      [srt{btype},ix{btype}] = sort(hipr{btype},'descend');
      vol{btype} = sort(ix{btype}(1:Ndots(btype)));
   end

   % x = 1:winL;
   % figure
   % plot(x(sort(ix(1:Ndots))),hipr(sort(ix(1:Ndots))),'o')

   for btype = 1:btypeN
      port1{btype}(vo,:) = (port{btype} + nrm(mean(f(Bnum{btype}(vol{btype}),:),1)))/2;%!!
   end
   subplot(121),plot(port1{1}(vo,:)),hold on
   subplot(122),plot(port1{2}(vo,:),'g'),hold on
end

%%
% Guessing
disp('Guessing')
win = 1:winL;

des = zeros(btypeN,btypeN,vom);
for vo = 1:vom
   disp(vo)
   for per = 1:perN

         for btype = 1:btypeN
            cor(btype,per) = f(per,win) * port1{btype}(vo,:)';
            cor(btype,per) = (cor(btype,per) +1)/2;
         end

      if ~any(per == Bnum{3})
         [~,ind] = max(cor(:,per));
         des(Bord(per),ind,vo) = des(Bord(per),ind,vo) + 1/Blen(Bord(per));
      end
   end
end

figure
k = 0;
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
      subplot(btypeN,btypeN,k),plot( permute(des(i,j,:),[1 3 2]) ,'.-'),axis tight,ylim([0 1])
  end
end
figure,plot(permute((des(1,1,:)+des(2,2,:))/2,[1 3 2]))

% figure
% plot(hipr)
% hold on,plot(cor(1,:),'r')
% hold on,plot(cor(2,:),'g')

