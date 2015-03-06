clc;
close all;
clear all;

load('D:\Dropbox\Signals\incartdb\I20\I20proc.mat')
in = val(1,:);

ann_num = 1:length(annot);
anniNAV{1} = ann_num(annot == 'N');  % Normal beat
anniNAV{2} = ann_num(annot == 'A');  % Atrial premature beat
anniNAV{3} = ann_num(annot == 'V');  % Premature ventricular contraction
ann{1} = mark(anniNAV{1});
ann{2} = mark(anniNAV{2});
ann{3} = mark(anniNAV{3});
ann_mark = 'NAV';
annNAV = ones(1,length(annot));
annNAV(ann_num(annot == 'A')) = 2*ones(1,length(ann{2}));
annNAV(ann_num(annot == 'V')) = 3*ones(1,length(ann{3}));

win = [-47 80];   % Borders of PQRST period
winL = win(2)-win(1)+1;
ki = 64;          % Interpolation coefficient
winXL = winL;%ki*(winL-1)+1;
qs = .08*.8;      % Quanting step (80% of P peak)
btypeN = 3;       % Beat types to examine


% Generating portraits
disp('Generating portraits')
perN = length(mark);  % Number of periods to use

f = zeros(perN, winL );
nf = zeros(perN, winL );
for per = 1:perN

   period = mark(per);
   window = period+win(1): period+win(2);
   f(per,:) = in(window);
   f(per,:) = f(per,:) - mean(f(per,:));

   nf(per,:) = nrm(f(per,:));

end

port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(anniNAV{btype},:),1));
end

dist = zeros(1,perN);
ndist = zeros(1,perN);
for per = 1:perN
   dist(per) = sum((port{annNAV(per)} - f(per,:)).^2).^0.5;
   ndist(per) = sum((port{annNAV(per)} - nf(per,:)).^2).^0.5;
end


%%
figure('color','white')
subplot(221)
plot(f(anniNAV{1}(1:40:end),:)',':k'),hold on
plot(port{1},'k','LineWidth',5),axis tight,ylim([-.6 1])
xlabel('a)')

nN = length(ann{1});
xval = .1:.1:2;

subplot(222)
[y,x] = hist(dist,xval);
stem(x,y/nN),grid
xlabel(sqrt(sum(dist.^2)/nN))

subplot(223)
plot(nf(anniNAV{1}(1:40:end),:)',':k'),hold on
plot(port{1},'k','LineWidth',5),axis tight,ylim([-.6 1])
xlabel('a)')

subplot(224)
[y,x] = hist(ndist,xval);
stem(x,y/nN),grid
xlabel(sqrt(sum(ndist.^2)/nN))



