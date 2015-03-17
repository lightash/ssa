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
btypeN = 3;       % Beat types to examine

% % Interpolating mean value
% f = zeros(length(annot),winL);
% m = zeros(1,length(annot));
% for per = 1:length(annot)
%    f(per,:) = in( mark(per)+win(1): mark(per)+win(2) );
%    m(per) = mean(f(per,:));
% end
% x = 1:size(in,2);
% mi = interp1(mark,m,x,'linear','extrap');
% in = in - mi;
% clear mi

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
xt = cell(1,perN);
ft = cell(1,perN);
for per = 1:perN

   period = mark(per);
   window = period+win(1): period+win(2);
   f(per,:) = in(window);
   f(per,:) = f(per,:) - mean(f(per,:));

   f(per,:) = nrm(f(per,:));

% %    x = 1:length(f(per,:));
% %    xi = linspace( 1, length(f(per,:)), ki*(length(f(per,:))-1)+1 );
% %    fi(per,:) = interp1(x,f(per,:),xi);
% % 
% %    [xt{per},ft{per}] = QRT(fi(per,:),qs,120,1);
end

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
% % 
% %    M(per) = mean(fqrt{annNAV(per)}(per,:));
% %    fqrt{annNAV(per)}(per,:) = nrm(fqrt{annNAV(per)}(per,:) - M(per));
% %    
% %    clear xs fs
% % end

port = cell(1,btypeN);
for btype = 1:btypeN
   port{btype} = nrm(mean(f(anniNAV{btype},:),1));
%    port{btype} = nrm( port{btype} - nrm(mean(f,1)) );
% %    port(btype,:) = nrm(mean(fqrt{btype},1));
end

%% Window width and position
disp('Window width and position')
tic

optwp = zeros(winL-1);
portwin = cell(winL-1);
for wid = 1:winL-1
   disp(wid)
   for pos = 1:winL-wid
      
      portwin{wid,pos} = pos:pos+wid;
      
      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      fg = f;
      for per = 1:perN
         
         for btype = 1:btypeN
            cor(btype,per) = nrm(fg(per,portwin{wid,pos}) - mean(fg(per,portwin{wid,pos})))*...
               nrm(port{btype}(portwin{wid,pos}) - mean(port{btype}(portwin{wid,pos})))';
         end
         
         [~,ind] = max(cor(:,per));
         des(annNAV(per),ind) = des(annNAV(per),ind) + 1/length(ann{annNAV(per)});
      end
      
      optwp(wid,pos) = des(1,1)+des(2,2)+des(3,3);
   end
end

toc
figure,surf(optwp),xlabel('Position'),ylabel('Width')

%% Informativity descent
disp('Informativity descent')
tic
cor = zeros(btypeN,perN);
des = zeros(btypeN);
indei = zeros(1,winL);
maxj = zeros(1,winL);
wini = cell(1,winL);
wini{1} = 1:winL;
% wini{1} = 1:2*winL;
fg = f;
% fg = [f f];
% for btype = 1:btypeN
%    port{btype} = [port{btype} port{btype}];
% end
for per = 1:perN
   for btype = 1:2%btypeN
      cor(btype,per) = nrm(fg(per,:) - mean(fg(per,:))) * ...
         nrm(port{btype} - mean(port{btype}))';
   end
   [~,ind] = max(cor(:,per));
   des(annNAV(per),ind) = des(annNAV(per),ind) + 1/length(ann{annNAV(per)});
end
% indei(1) = (des(1,1)+des(2,2)+des(3,3))/3;
indei(1) = (des(1,1)+des(2,2))/2;

for i = 2:winL
% for i = 2:2*winL
   disp(i)

   indej = zeros(1,length(wini{i-1}));
   for j = 1:length(wini{i-1})

      winj = [wini{i-1}(1:j-1) wini{i-1}(j+1:end)];

      cor = zeros(btypeN,perN);
      des = zeros(btypeN);
      for per = 1:perN
         
         for btype = 1:2%btypeN
            cor(btype,per) = nrm(fg(per,winj) - mean(fg(per,winj))) * ...
               nrm(port{btype}(winj) - mean(port{btype}(winj)))';
         end
         
         [~,ind] = max(cor(:,per));
         des(annNAV(per),ind) = des(annNAV(per),ind) + 1/length(ann{annNAV(per)});
         
      end
%       indej(j) = (des(1,1)+des(2,2)+des(3,3))/3;
      indej(j) = (des(1,1)+des(2,2))/2;
      
   end
   
   [indei(i),maxj(i)] = max(indej);

   wini{i} = [wini{i-1}(1:maxj(i)-1) wini{i-1}(maxj(i)+1:end)];
      
end

toc
save('indei_NA','indei','wini','maxj')
x = 1:winL;
% x = 1:2*winL;
[~,maxi] = max(indei);
figure,plot(x,port{1},'-',x,port{2},'-g',x,port{3},'-r')
hold on,plot(wini{maxi-1},port{1}(wini{maxi-1}),'.',wini{maxi-1},port{2}(wini{maxi-1}),'.g',...
   wini{maxi-1},port{3}(wini{maxi-1}),'.r'),grid,axis tight
figure,subplot(211),plot(indei,'.-'),axis tight,subplot(212),plot(maxj,'.-'),axis tight
%% Guessing
disp('Guessing')
% [wid,pos] = find( optwp == max(max(optwp)) );
% winwp = pos:pos+wid;
% % perN = length(mark);
% % for per = 1:perN
% %    
% %    period = mark(per);
% %    window = period+win(1): period+win(2);
% %    f = in(window);
% %    f = f - mean(f);
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

% partN = winL/4;
% partS = floor(winXL/partN);
% cor = cell(1,partN);
% des = cell(1,partN);
% for parti = 1:partN
%    des{parti} = zeros(btypeN);
% end
des = zeros(btypeN);
% % fg = zeros(perN,winXL);
fg = f;
for per = 1:perN
   disp(per)
% %    x = 1:length(xi);
% %    fg(per,:) = interp1(xt{per},ft{per},x,'pchip','extrap');
% %    fg(per,:) = fg(per,:)-mean(fg(per,:));
   
%    for parti = 1:partN%[1:6 8:partN]%[1:24 26:partN]
      
%       if parti > 7
%          window = [1:6*partS 8:parti*partS];
%       else
%          window = 1:parti*partS;%winXL;(parti-1)*partS+
%       end
%       fg(per,:) = fg(per,:) + max(fg(per,:))*randn(1,winL)/12;
      sig(per,:) = (nrm(fg(per,:) - mean(fg(per,:))));% - nrm(mean(f,1)));
      for btype = 1:btypeN
         cor(btype,per) = sig(per,:) * port{btype}';
%          cor(btype,per) = nrm(fg(per,winwp) - mean(fg(per,winwp))) * nrm(port{btype}(winwp) - mean(port{btype}(winwp)))';%mean(fg(per,:))mean(port{btype})
%          a = nrm(port{btype}(1:2:end) - mean(port{btype}(1:2:end)));
%          b = nrm(port{btype}(2:2:end) - mean(port{btype}(2:2:end)));
%          c = nrm(fg(per,1:2:end) - mean(fg(per,1:2:end)));
%          d = nrm(fg(per,2:2:end) - mean(fg(per,2:2:end)));
%          cor(btype,per) = sum( a.^2+b.^2 .* a.*c + b.*d );
         cor(btype,per) = (cor(btype,per) +1)/2;
      end
      
      [~,ind] = max(cor(:,per));
%       if ~any(per == anniNAV{1})
      des(annNAV(per),ind) = des(annNAV(per),ind) + 1/length(ann{annNAV(per)});
%       end
%    end
end

%%
figure
k = 0;
% des_m = zeros(btypeN*btypeN,partN);
% ddes = zeros(btypeN*btypeN,partN-1);
for i = 1:btypeN
   for j = 1:btypeN
      k = k+1;
      
%       for parti = 1:partN
%          des_m(k,parti) = des{parti}(i,j);
%       end
%       ddes(k,:) = diff(des_m(k,:));
      
      subplot(btypeN,btypeN,k),stem(des(i,j),'.-'),axis([0 2 0 1])%,grid%+btypeN*(i-1)
      xlabel(des(i,j))
%       subplot(6,btypeN,k+btypeN*i),plot(ddes(k,:),'.-'),grid,axis tight
  end
end
title((des(1,1)+des(2,2)+des(3,3))/3)
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






