clc;
% close all;
clear all;

load('D:\Dropbox\Signals\incartdb\I20\I20proc.mat')
in = val(1,:);
% load('D:\Dropbox\Signals\incartdb\I21\I21proc.mat')
% in20 = val(1,:);

ann_num = 1:length(annot);
ann{1} = mark(ann_num(annot == 'N'));  % Normal beat
ann{2} = mark(ann_num(annot == 'A'));  % Atrial premature beat
ann{3} = mark(ann_num(annot == 'V'));  % Premature ventricular contraction
ann_mark = 'NAV';

win = [-47 80];   % Borders of PQRST period
winL = win(2)-win(1)+1;
ki = 64;          % Interpolation coefficient
qs = .08*.8;      % Quanting step

% Interpolating mean value
f = zeros(length(annot),winL);
m = zeros(1,length(annot));
for per = 1:length(annot)
   f(per,:) = in( mark(per)+win(1): mark(per)+win(2) );
   m(per) = mean(f(per,:));
end
x = 1:size(in,2);
mi = interp1(mark,m,x,'linear','extrap');
in = in - mi;
clear mi

% % Filtering
% h = fdesign.lowpass('N,Fc', 20, 40/Fd*2);
% d = design(h);
% in_f = filtfilt(d.Numerator,1,in);
% % freqz(d)
% figure,plot(in,'r'),hold on,plot(in_f,'b'),ylim([-.5 .5])
% figure,plot(in-in_f)

% Generating portraits
for btype = 1:3  % Beat type
   disp(btype)
   perN = length(ann{btype});  % Number of periods to use
   
   f = zeros(perN, win(2)+1-win(1) );
   fi = zeros(perN, ki*(win(2)-win(1)) );
   fqrt = zeros(perN, ki*(win(2)-win(1))+1 );%{btype}
   xqrt{btype} = [];
   xt_all = cell(1,perN);
   ft_all = cell(1,perN);
   for per = 1:perN

      period = ann{btype}(per);
      window = period+win(1): period+win(2);
      f(per,:) = in(window);
      f(per,:) = f(per,:) - mean(f(per,:));

%       % Filtering
%       filtsize = 4;
%       f(per,:) = filtfilt(ones(1,filtsize)/filtsize,1,f(per,:));
   %    h = fdesign.highpass('N,Fc', 50, 3/Fd*2);
   %    d = design(h);
   %    f = filtfilt(d.Numerator,1,f);
   %    freqz(d)

      f(per,:) = nrm(f(per,:));

      x = 1:length(f(per,:));
      xi = linspace( 1, length(f(per,:)), ki*(length(f(per,:))-1)+1 );
      fi = interp1(x,f(per,:),xi);

      [xt{btype,per},ft{btype,per}] = QRT(fi,qs,120);
      
   end
end
   
% Interpolating mean value
m = zeros(1,length(mark));
for per = 1:length(mark)
   m(per) = mean(ft{btype,per});
end
x = 1:size(in,2);
mi(1,:) = interp1([mark(1) mark(2)],[m(1) m(2)],...
   ki*(mark(1)+win(1)-1)+1:ki*(mark(1)+win(2)-1)+1,'linear','extrap');
for per = 2:perN-1
   disp(per)
   left = interp1([mark(per-1) mark(per)],[m(per-1) m(per)],ki*(mark(per)+win(1)-1)+1:ki*(mark(per)-1)+1);
   right = interp1([mark(per) mark(per+1)],[m(per) m(per+1)],ki*(mark(per)-1)+1:ki*(mark(per)+win(2)-1)+1);
   mi(per,:) = [left right(2:end)];
end
mi(perN,:) = interp1([mark(end-1) mark(end)],[m(end-1) m(end)],...
   ki*(mark(end)+win(1)-1)+1:ki*(mark(end)+win(2)-1)+1,'linear','extrap');
figure,plot(m,'.-'),hold on,plot(mi,'r')

for btype = 1:3
   for per = 1:perN
%       period = (ann{btype}(per)-1)*ki+1;
%       window = period+ki*win(1) : period+ki*win(2);
      ft{per} = ft{per} - mi(window);
   end
   
   for per = 1:perN
      
      x = 1:length(xi);
      fqrt(per,:) = interp1(xt,ft,x);
      
      fqrt(per,:) = nrm(fqrt(per,:) - mean(fqrt(per,:)));
      
%       [xs,fs] = stairs(xt,[ft(2:end) ft(end)]);
%       fqrt(per,1) = 0;
%       for i = 2:length(xs)
%          fqrt(per,xs(i-1)+1:xs(i)) = fs(i);
%       end
%       fqrt(per,end) = 0;
% %       xt_all{per} = xt;
% %       ft_all{per} = ft;
% %       xqrt{btype} = [xqrt{btype} xt];
      
      clear fi xt ft xs fs
   end
   
%    xqrt{btype} = unique(sort(xqrt{btype}));
%    fqrt{btype} = zeros(perN,length(xqrt{btype}));
%    for per = 1:perN
%       fqrt{btype}(per,:) = interp1(xt_all{per},ft_all{per},xqrt{btype});
%    end
   
%    figure
%    subplot(211),plot(f','-b'),grid,axis tight,hold on
%    subplot(212),stem(x,fqrt','.b'),axis tight,grid on,hold on
   
%    port(btype,:) = nrm(mean(f,1));
   port(btype,:) = nrm(mean(fqrt,1));%{btype}
end
% figure,plot(port')
% load(['D:\Dropbox\Signals\incartdb\I20\I20port.mat'])
%% Guessing
perN = length(mark)-1;%size(val,2)-winL;
cor = zeros(3,perN);
des = zeros(3);
for per = 1:perN
   
   period = mark(per);
   window = period+win(1): period+win(2);
   f = in(window);
   f = nrm(f - mean(f));
   
   x = 1:length(f);
   xi = linspace(1, length(f), ki*length(f));
   fi = interp1(x,f,xi);
   
   [xt,ft{per}] = QRT(fi,qs,120);
   end
   
   % Interpolating mean value
   m = zeros(1,length(annot));
   for per = 1:length(annot)
      m(per) = mean(ft{per});
   end
   x = 1:size(in,2);
   mi = interp1(mark,m,x,'linear','extrap');

   for per = 1:perN
      period = mark(per);
      window = period+win(1): period+win(2);
      ft{per} = ft{per} - mi(window);
   end
   
   for per = 1:perN
   ft{per} = nrm(ft{per}-mean(ft{per}));
   
   x = 1:length(xi);
   fg = interp1(xt,ft{per},x);
%    [xs,fs] = stairs(xt,[ft(2:end) ft(end)]);
%    fg(1) = 0;
%    for i = 2:length(xs)
%       fg(xs(i-1)+1:xs(i)) = fs(i);
%    end
%    fg(end) = 0;
   
   clear fi xt ft xs fs
   
   for btype = 1:3
      cor(btype,per) = fg * port(btype,:)';
   end
   
   [~,ind] = max(cor(:,per));
   if annot(per) == 'N'
      des(1,ind) = des(1,ind)+1;
   elseif annot(per) == 'A'
      des(2,ind) = des(2,ind)+1;
   elseif annot(per) == 'V'
      des(3,ind) = des(3,ind)+1;
   end
   
end
%%
figure
subplot(311),stem(des(1,:)/length(ann{1})),axis([0 4 0 1])
xlabel(num2str(des(1,:)/length(ann{1})))
subplot(312),stem(des(2,:)/length(ann{2})),axis([0 4 0 1])
xlabel(num2str(des(2,:)/length(ann{2})))
subplot(313),stem(des(3,:)/length(ann{3})),axis([0 4 0 1])
xlabel(num2str(des(3,:)/length(ann{3})))

%%
% figure
% stem(ann{1},ones(1,length(ann{1})),'.:r'),hold on
% stem(ann{2},.9*ones(1,length(ann{2})),'o:r'),hold on
% stem(ann{3},1.1*ones(1,length(ann{3})),'x:r'),hold on
% plot(mark(1:end-1),cor,'-'),ylim([.2 1.2])%-win(1):perN-win(1),
% plot(mark(1:end-1),des1,'.-'),hold on,plot(mark(1:end-1),des2,'.-r'),grid









