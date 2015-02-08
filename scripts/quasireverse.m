clc;
% close all;
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
btypeN = 3;

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
disp('Generating portraits')
perN = length(mark);  % Number of periods to use

f = zeros(perN, winL );
fi = zeros(perN, winXL );
% xqrt = [];
xt = cell(1,perN);
ft = cell(1,perN);
for per = 1:perN

   period = mark(per);
   window = period+win(1): period+win(2);
   f(per,:) = in(window);
   f(per,:) = f(per,:) - mean(f(per,:));

%       % Filtering
%       filtsize = 4;
%       f(per,:) = filtfilt(ones(1,filtsize)/filtsize,1,f(per,:));
%    h = fdesign.highpass('N,Fc', 50, btypeN/Fd*2);
%    d = design(h);
%    f = filtfilt(d.Numerator,1,f);
%    freqz(d)

   f(per,:) = nrm(f(per,:));

% %    x = 1:length(f(per,:));
% %    xi = linspace( 1, length(f(per,:)), ki*(length(f(per,:))-1)+1 );
% %    fi(per,:) = interp1(x,f(per,:),xi);
% % 
% %    [xt{per},ft{per}] = QRT(fi(per,:),qs,120,1);
end

% for btype = 1:btypeN
%    disp(btype)
%    port_s{btype} = fi(anniNAV{btype},:);
%    Me = mean(port_s{btype},2);
%    for per = 1:size(port_s{btype},1)
%       port_s{btype}(per,:) = port_s{btype}(per,:) - Me(per);
%    end
%    minmax(btype,:) = [mean(min(port_s{btype},[],1)) mean(max(port_s{btype},[],1))];
% end

% % % Interpolating mean value in portraits
% % m = zeros(1,length(mark));
% % xt_all = [];
% % for per = 1:length(mark)
% %    m(per) = mean(interp1(xt{per},ft{per},1:winXL,'linear','extrap'));
% %    xt_all = [xt_all xt{per}+ki*(mark(per)+win(1)-1)+1];
% % end
% % mi = interp1(ki*(mark-1)+1,m,xt_all,'linear','extrap');
% % for per = 1:perN
% %    ft{per} = ft{per} - mi(1:length(xt{per}));
% %    mi = mi(length(xt{per})+1:end);
% % end
% % 
% % fqrt = cell(1,btypeN);
% % for btype = 1:btypeN
% %    fqrt{btype} = zeros(perN,winXL);
% % end
% % M = zeros(1,perN);
% % for per = 1:perN
% %    disp(per)
% %    
% %    x = 1:length(xi);
% %    fqrt{annNAV(per)}(per,:) = interp1(xt{per},ft{per},x,'pchip','extrap');
% % %    [xs,fs] = stairs(xt{per},[ft{per}(2:end) ft{per}(end)]);
% % %    fqrt{annNAV(per)}(per,1) = 0;
% % %    for i = 2:length(xs)
% % %       fqrt{annNAV(per)}(per,xs(i-1)+1:xs(i)) = fs(i);
% % %    end
% % %    fqrt{annNAV(per)}(per,end) = 0;
% % % %       xt_all{per} = xt;
% % % %       ft_all{per} = ft;
% % % %       xqrt{btype} = [xqrt{btype} xt];
% % 
% %    M(per) = mean(fqrt{annNAV(per)}(per,:));
% %    fqrt{annNAV(per)}(per,:) = nrm(fqrt{annNAV(per)}(per,:) - M(per));
% %    
% %    clear xs fs
% % end

%    xqrt{btype} = unique(sort(xqrt{btype}));
%    fqrt{btype} = zeros(perN,length(xqrt{btype}));
%    for per = 1:perN
%       fqrt{btype}(per,:) = interp1(xt_all{per},ft_all{per},xqrt{btype});
%    end

port = zeros(btypeN,winL);
for btype = 1:btypeN
   port(btype,:) = nrm(mean(f(anniNAV{btype},:),1));
% %    port(btype,:) = nrm(mean(fqrt{btype},1));
end

% load(['D:\Dropbox\Signals\incartdb\I20\I20port.mat'])
%% Guessing
disp('Guessing')
% % perN = length(mark);
% % for per = 1:perN
% %    
% %    period = mark(per);
% %    window = period+win(1): period+win(2);
% %    f = in(window);
% %    f = f - mean(f);
% %    
% % %    f = f*min(abs( [minmax(1)/min(f) minmax(2)/max(f)] ));
% %    
% %    f = nrm(f);
% %    x = 1:length(f);
% %    xi = linspace(1, length(f), ki*(length(f)-1)+1);
% %    fi = interp1(x,f,xi);
% % 
% %    [xt{per},ft{per}] = QRT(fi,qs,120,1);
% % end
% % 
% % % Interpolating mean value for guessing
% % m = zeros(1,perN);
% % for per = 1:perN
% %    m(per) = mean(interp1(xt{per},ft{per},1:winXL,'linear','extrap'));
% %    xt_all = [xt_all xt{per}+ki*(mark(per)+win(1)-1)+1];
% % end
% % mi = interp1(ki*(mark-1)+1,m,xt_all,'linear','extrap');
% % for per = 1:perN
% %    ft{per} = ft{per} - mi(1:length(xt{per}));
% %    mi = mi(length(xt{per})+1:end);
% % end

partN = winL/4;
partS = floor(winXL/partN);
cor = cell(1,partN);
des = cell(1,partN);
for parti = 1:partN
   des{parti} = zeros(btypeN);
end
% % fg = zeros(perN,winXL);
fg = f;
for per = 1:perN
   disp(per)
% %    x = 1:length(xi);
% %    fg(per,:) = interp1(xt{per},ft{per},x,'pchip','extrap');
%    [xs,fs] = stairs(xt{per},[ft{per}(2:end) ft{per}(end)]);
%    fg(1) = 0;
%    for i = 2:length(xs)
%       fg(xs(i-1)+1:xs(i)) = fs(i);
%    end
%    fg(end) = 0;
   
% %    fg(per,:) = fg(per,:)-mean(fg(per,:));
%    clear xs fs
   
   for parti = 1:partN%[1:6 8:partN]%[1:24 26:partN]
      
%       if parti > 7
%          window = [1:6*partS 8:parti*partS];
%       else
         window = 1:parti*partS;%winXL;(parti-1)*partS+
%       end

      for btype = 1:btypeN
         cor{parti}(btype,per) = nrm(fg(per,window) - mean(fg(per,window))) * nrm(port(btype,window) - mean(port(btype,window)))';
      end
      
      [~,ind] = max(cor{parti}(:,per));
      if ~any(per == anniNAV{1})
         des{parti}(annNAV(per),ind) = des{parti}(annNAV(per),ind)+...
            1/length(ann{annNAV(per)});
      end
   end
end

figure
k = 0;
des_m = zeros(btypeN*btypeN,partN);
ddes = zeros(btypeN*btypeN,partN-1);
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
      for parti = 1:partN
         des_m(k,parti) = des{parti}(i,j);
      end
      ddes(k,:) = diff(des_m(k,:));
      
      subplot(6,btypeN,k+btypeN*(i-1)),plot(des_m(k,:),'.-'),axis([1 partN+1 0 1]),grid
      subplot(6,btypeN,k+btypeN*i),plot(ddes(k,:),'.-'),grid,axis tight
  end
end
% %%
% % Sorting
% btype = 1;
% [~,Srt_ind] = sort(ddes(btype,:),'descend');
% 
% Srt = zeros(perN,winXL);
% Srt_cor = cell(1,partN);
% Srt_des = zeros(btypeN,btypeN,partN);
% Srt_des = num2cell(Srt_des,[1 2]);
% for per = 1:perN
%    for ai = 1:length(Srt_ind)
%       Srt(per,partS*(ai-1)+1:partS*ai) = fg(per,partS*(Srt_ind(ai)-1)+1:partS*Srt_ind(ai));
%    end
%    
%    for parti = 1:partN
%       window = 1:parti*partS;
%       
%       for btype = 1:btypeN
%          Srt_cor{parti}(btype,per) = nrm(Srt(per,window)) * nrm(port(btype,window))';
%       end
%       
%       [~,ind] = max(Srt_cor{parti}(:,per));
%       Srt_des{parti}(annNAV(per),ind) = Srt_des{parti}(annNAV(per),ind)+...
%          1/length(ann{annNAV(per)});
%       
%    end
% end
% 
% figure
% k = 0;
% Srt_des_m = zeros(btypeN*btypeN,partN);
% Srt_ddes = zeros(btypeN*btypeN,partN-1);
% for i = 1:btypeN
%    for j = 1:btypeN
%       k = k+1;
%       
%       for parti = 1:partN
%          Srt_des_m(k,parti) = Srt_des{parti}(i,j);
%       end
%       Srt_ddes(k,:) = diff(Srt_des_m(k,:));
%       
%       subplot(6,btypeN,k+btypeN*(i-1)),plot(Srt_des_m(k,:),'.-'),axis([1 partN+1 0 1])
%       subplot(6,btypeN,k+btypeN*i),plot(Srt_ddes(k,:),'.-'),axis tight
%   end
% end
% 






