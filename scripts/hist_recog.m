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

btypeN = 2;

% Generating histograms
Nbars = 16;
for i = 1:winL
   [~,scale(i,:)] = hist( [f(Bnum{1},i);f(Bnum{2},i)] ,Nbars);
   
   for btype = 1:btypeN
      Mh{btype}(i) = mean(f(Bnum{btype},i));
      H{btype}(:,i) = hist( f(Bnum{btype},i) ,scale(i,:))/Blen(btype)';
      % Распределить наименьший столбец по нулям
      [minh,indh] = min(H{btype}( H{btype}(:,i)>0, i));
      mzh = minh/( length( H{btype}( H{btype}(:,i)==0, i)) +1);
      H{btype}( H{btype}(:,i)==0 | H{btype}(:,i)==minh, i) = mzh;
      H{btype}(:,i) = H{btype}(:,i)/sum(H{btype}(:,i));
   end
end
%%
% Border between p00 and p11

% border = zeros(1,winL);
for i = 1:winL
   m1 = Mh{1}(i);
   m2 = Mh{2}(i);
   if m1 < m2
      leftN(i) = 1;
   else
      leftN(i) = 0;
   end
   
   Nbars = 256;
   scali(i,:) = linspace(scale(i,1),scale(i,end),Nbars);
   for btype = 1:btypeN
      Hi{btype}(:,i) = interp1(scale(i,:),H{btype}(:,i),scali(i,:));
      Hi{btype}(:,i) = Hi{btype}(:,i)/sum(Hi{btype}(:,i));
   end
   
   for bari = 1:Nbars
      if leftN(i)
         brd(bari) = sum(Hi{1}( 1:bari ,i)) - sum(Hi{2}( bari+1:end ,i));
      else
         brd(bari) = sum(Hi{2}( 1:bari ,i)) - sum(Hi{1}( bari+1:end ,i));
      end
   end
   [~,border(i)] = min(abs(brd));
   
   if leftN(i)
      p00(i) = sum(Hi{1}(1:border(i),i));
      p11(i) = sum(Hi{2}(border(i)+1:end,i));
   else
      p00(i) = sum(Hi{1}(border(i)+1:end,i));
      p11(i) = sum(Hi{2}(1:border(i),i));
   end
end
% %%
% for i = 1:10:winL
%    figure(i)
%    plot(scali(i,:),Hi{1}(:,i),'.-'),hold on
%    plot(scali(i,:),Hi{2}(:,i),'.-g'),hold on
%    plot([scali(i,border(i)) scali(i,border(i))]+.5*mean(diff(scali(i,:))),[0 max([Hi{1}(:,i);Hi{2}(:,i)])],'k'),grid
%    title(num2str([sum(Hi{1}(1:border(i),i)) sum(Hi{2}(border(i)+1:end,i)) leftN(i)]))
% end
%%
figure
plot(diff([p00;p11]),'.-'),grid,axis([1 winL -1 1])














