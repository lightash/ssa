clc;close all;clear all;

load('D:\Dropbox\MATLAB\Био. комп. сис\signal\Ms.mat')

N = 64;
n = 16;
A = Ms(1:N,1:n:n*N);

noi = 1e-6:.2:1e-6 + 1;  % Noise vector
stn = 5%00;               % Statistics size
ord = -1:.2:1;           % Histogram bars

[a, b, e] = ABAxTriDecomp(A, 0, noi, stn, ord);%, 'wrkspc_13-12-08_x1.mat');

A2 = zeros(2*N);
for i = 1:N-1
   A2(i,:) = [A(i,:) A(i+1,:)];
end
A2(N,:) = [A(N,:) A(1,:)];
j = 0;
for i = N+1:2*N-2
   j = j+1;
   A2(i,:) = [A(j,:) A(j+2,:)];
end
A2(2*N-1,:) = [A(N-3,:) A(N-1,:)];
A2(2*N,:) = [A(N-2,:) A(N,:)];

[a2, b2] = ABAxTriDecomp(A2, noi, stn, ord);%, 'wrkspc_13-12-08_x2.mat');

%%

figure
plot(noi,a.am_cor_a,'.k',noi,a.am_cor_tri,'+:g',noi,a.am_cor_ax,'x-r',noi,a.am_cor_b,'*:b')
axis tight
% ...for double components
   figure
   plot(noi,a2.am_cor_a,'.k',noi,a2.am_cor_tri,'+:g',noi,a2.am_cor_ax,'x-r',noi,a2.am_cor_b,'*:b')
   axis tight

%%

% clc;
% close all;
% clear all;
% 
% load('wrkspc_13-12-08_x1.mat')
% load('wrkspc_13-12-08_x2.mat')



[a, b, c] = ZeroAltErrProbs(a, b, noi, ord);
[a2, b2, c2] = ZeroAltErrProbs(a2, b2, noi, ord);

% n = noiN;
% m = 1;
% % for n = 1:noiN
%    d(1,:) = a.hist_ax(n,m,:);f(1,:) = b.hist_ax(n,m,:);
%    figure,subplot(211),plot(ord,d,'xb'),hold on,plot(ord,f,'+g')
% % end


% for i = 1:N
%    figure,plot(noi,ord(c.T_ax(:,1)),'.-'),axis tight
% end


% close all
% n = noiN;
% m = 21;
% for n=1:noiN
%    var1(1,:) = c.P1_ax(n,m,:);var2(1,:) = c.P2_ax(n,m,:);
%    figure,plot(var1,'x-b'),hold on,plot(var2,'+-g'),axis tight
% end


figure,plot(noi,c.p00_a(:,1),'-r',noi,c.p01_a(:,1),':b',noi,c.p10_a(:,1),':r',noi,c.p11_a(:,1),'-b'),axis tight
legend('00','01','10','11')

figure,plot(noi,c.p00_ax(:,1),'-r',noi,c.p01_ax(:,1),':b',noi,c.p10_ax(:,1),':r',noi,c.p11_ax(:,1),'-b'),axis tight
legend('00','01','10','11')


figure,plot(noi,c2.p00_a(:,1),'-r',noi,c2.p01_a(:,1),':b',noi,c2.p10_a(:,1),':r',noi,c2.p11_a(:,1),'-b'),axis tight
legend('00','01','10','11')

figure,plot(noi,c2.p00_ax(:,1),'-r',noi,c2.p01_ax(:,1),':b',noi,c2.p10_ax(:,1),':r',noi,c2.p11_ax(:,1),'-b'),axis tight
legend('00','01','10','11')






















