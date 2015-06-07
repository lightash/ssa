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

Nbins = 16;
clear hyp scale H
for i = 1:winL
   [hyp(i,:),scale(i,:)] = hist( [f(Bnum{1},i);f(Bnum{2},i)] ,Nbins);
   
   for btype = 1:btypeN
      H{btype}(:,i) = hist( f(Bnum{btype},i) ,scale(i,:))/Blen(btype)';
      % Распределить наименьший столбец по нулям
      [minh,indh] = min(H{btype}( H{btype}(:,i)>0, i));
      mzh = minh/( length( H{btype}( H{btype}(:,i)==0, i)) +1);
      H{btype}( H{btype}(:,i)==0 | H{btype}(:,i)==minh, i) = mzh;
   end
end
I(1,:) = DK( H{1}, H{2}, Blen(1), Blen(2) );
I(2,:) = AlphaZ( H{1}, H{2} );

% Energy & dispersion
I(3,:) = ED(f(Bnum{1},:),f(Bnum{2},:),1);
[~,I(4,:)] = ED(f(Bnum{1},:),f(Bnum{2},:),1);
I(5,:) = ED(f(Bnum{1},:),f(Bnum{2},:),2);
[~,I(6,:)] = ED(f(Bnum{1},:),f(Bnum{2},:),2);
I(7,:) = ED(f(Bnum{1},:),f(Bnum{2},:),3);
[~,I(8,:)] = ED(f(Bnum{1},:),f(Bnum{2},:),3);

I(9,:) = HB(f(Bnum{1},:),f(Bnum{2},:));

load('indei_NA.mat','inform')
I(10,:) = inform;

for form = 1:size(I,1)
   disp(form)
   Ikl = I(form,:);

   x = 1:winL;

   for Ndots = 1:winL
      % Ndots = 100;%length(wini{maxi-1});
      [srt,ix] = sort(Ikl,'descend');

      % Guessing
      win = sort(ix(1:Ndots));

      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      for per = 1:perN
         sig(per,1:Ndots) = (nrm(f(per,win) - mean(f(per,win))));
         for btype = 1:btypeN
            cor(btype,per) = sig(per,1:Ndots) * port{btype}(win)';
            cor(btype,per) = (cor(btype,per) +1)/2;
         end

         [~,ind] = max(cor(:,per));
         if ~any(per == Bnum{3})
            des(Bord(per),ind) = des(Bord(per),ind) + 1/length(Bpos{Bord(per)});
         end
      end
      des_m(form,Ndots) = (des(1,1)+des(2,2))/2;
   end

% figure
% subplot(2,2,form)
% plot(Ikl),grid,axis tight
% xlabel(num2str([des(1,1) des(2,2) (des(1,1)+des(2,2))/2]))
% hold on,plot(x(sort(ix(1:Ndots))),Ikl(sort(ix(1:Ndots))),'.')

subplot(5,2,form)
plot(des_m(form,:),'.-'),grid,axis([1 winL .5 1])

end










