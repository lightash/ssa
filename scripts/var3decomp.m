clc;
close all;
clear all;

load('g:\Dropbox\MATLAB\Био. комп. сис\signal\Ms.mat')
N = 64;%32;
n = 16;%32;
si = transform(.1*sin(2*pi*50*[0:1/N^2:1]),'matrix');
A = Ms(1:N,1:n:n*N)+0*randn(N) + si;

[a2,b] = OSR(A);

acp = mean(A,1);
l = acp - b;

ai = a2 - transform(l,'matrix_repeat');

for i=1:N
    aipl(i,:) = ( ai(i,:) * (l/sqrt(l*l'))' ) * l/sqrt(l*l');
    aiol(i,:) = ai(i,:) - aipl(i,:);
end

figure('Color','w')
subplot(711),plot(A'),axis tight,grid on,ylabel('A')
subplot(712),plot(b),grid on,axis tight,ylabel('b')
subplot(713),plot(a2'),grid on,axis tight,ylabel('a2')
subplot(714),plot(l),grid on,axis tight,ylabel('l')
subplot(715),plot(ai'),grid on,axis tight,ylabel('ai')
subplot(716),plot(aipl'),grid on,axis tight,ylabel('aipl')
subplot(717),plot(aiol'),grid on,axis tight,ylabel('aiol')

figure,plot(transform(A,'vector'),'g')
hold on,plot(transform(l,'vector_repeat'),'.-r')
hold on,plot(transform(Ms(1:N,1:n:n*N),'vector'),'b')

% for i = 1:N
%     loai(i) = l*aiol(i,:)';
% end
% disp(max(abs(loai)))
% 
% for i = 1:N
%     acpoai(i) = acp*aiol(i,:)';
% end
% disp(max(abs(acpoai)))
% 
% for i = 1:N
%     for j = 1:N
%         aiop(i,j) = aipl(i,:)*aiol(j,:)';
%     end
% end
% disp(max(max(abs(aiop))))




