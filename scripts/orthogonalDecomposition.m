clc;
% close all;
clear all;

N1 = 16;
x = 215:5:290;
S1 = [1.47 1.08 .77 .53  .37 .32 .31 .33  .35 .37 .4  .35  .3 .26 .2 .16];
S2 = [1.2 1.03 .88 .78  .72 .71 .8 .93  1.03 1  .9 .76  .68 .47 .3 .21];

% P11 = zeros(1,11);
% for sigma = 0:0.1:1
%     right = 0;
%     for i = 1:1:100
% 
%         noise = sigma*randn(1,N1);

ps1 = sqrt(S1*S1');     % length of S1
ps2 = sqrt(S2*S2');     %           S2

e1 = S1/ps1;            % orth of S1
e2 = S2/ps2;            %         S2

ps12 = S1*e2';          % length of projection of S1 on e2
ps21 = S2*e1';          %                         S2    e1

S12 = ps12*e2;          % projection of S1 on e2
S21 = ps21*e1;          %               S2    e1

S11 = S1 - S12;         % exclusive component of S1
S22 = S2 - S21;         %                        S2

ps11 = sqrt(S11*S11');  % length of S11
ps22 = sqrt(S22*S22');  %           S22

% figure,plot(x,S1,'r',x,S2,'b',x,S11,'--r',x,S22,'--b',x,S12,'-.r',x,S21,'-.b')
% legend('S1','S2','S11','S22','S12','S21'),grid on,axis tight

xi = linspace(x(1),x(end),160);
S1i = spline(x,S1,xi);
S2i = spline(x,S2,xi);
S11i = spline(x,S11,xi);
S22i = spline(x,S22,xi);
S12i = spline(x,S12,xi);
S21i = spline(x,S21,xi);

figure('Color','w')
plot(xi,S1i,'k',xi,S2i,'k',xi,S11i,'-.k',xi,S22i,'-.k',xi,S12i,':k','LineWidth',1.5),hold on,plot(xi,S21i,'k')
legend('S_1','S_2','S_1_1','S_2_2','S_1_2','S_2_1'),grid on,axis tight
hold on,plot(x,S1,'.k',x,S2,'.k',x,S11,'.k',x,S22,'.k',x,S12,'.k',x,S21,'.k','MarkerSize',15)
xlabel('Wavelength \lambda, nm','FontName','Times New Roman','FontSize',10)
ylabel('The extinction coefficient  \aleph, cm^-^1','FontName','Times New Roman','FontSize',10)


% pe(1,1) = ps11/ps1;
% pe(1,2) = ps12/ps1;
% pe(2,1) = ps21/ps2;
% pe(2,2) = ps22/ps2;
% 
% dS = S2 - S1;
% beta = 0.2;
% x = S1 + beta*dS;
% teta = 0.8;
% x = teta*x + noise;
% 
% px = sqrt(x*x');
% ex = x/px;
% 
% pex = zeros(2);
% pex(1,2) = sqrt(ex*e2');
% pex(1,1) = sqrt( (ex - pex(1,2)*e2)*(ex - pex(1,2)*e2)' );
% pex(2,1) = sqrt(ex*e1');
% pex(2,2) = sqrt( (ex - pex(2,1)*e1)*(ex - pex(2,1)*e1)' );
% 
% r = 0.5*(1 - abs(pe - pex)/2);
% 
% rr1 = sum(r(1,:));
% rr2 = sum(r(2,:));
% 
% if rr1 >= rr2
% %             right = right+1;
%     disp('S1')
% else
%     disp('S2')
% end
% 
%     end
%     P11(sigma*10+1) = right/100;
% end
% 
% figure,plot(P11)
% S1_1 = S1 - S1_11;
% S2_1 = S2 - S2_11;
% 
% disp(S1_1*S1_11')
% disp(S2_1*S2_11')
% 
% figure,plot(1:N1,S1,1:N1,S1_1,'g',1:N1,S1_11,'r'),axis([1 16 -1 2]),grid
% figure,plot(1:N1,S2,1:N1,S2_1,'g',1:N1,S2_11,'r'),axis([1 16 -1 2]),grid
