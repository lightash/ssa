clc;
close all;
clear all;

load('Signals\Msasr')
N1 = 16;
s = 1;%503;
n = 50;%1;
a1 = Msasr(1:N1,s:n:s+n*N1-1);

[a_1, a_coll, a_ort, S] = preprocess(a1, 0);

a_c = transform(a_coll,'vector');
a_o = transform(a_ort,'vector');
aa = transform(a_1,'vector');
ss = transform(S,'vector_repeat');

figure
subplot(311),plot(aa,'k'),axis tight,title('Исходная последовательность a'),grid
subplot(312),plot(1:N1*N1,a_c,'--k',1:N1*N1,a_o,'k'),axis tight,title('Колл., орт. составляющие'),grid
subplot(313),plot(ss,'k'),axis tight,title('Среднее'),grid

S = ss-mean(ss);%aa;% + .5*randn(1,length(aa));
Et = ss-mean(ss);disp(mean(ss))
N = N1*N1;
NN = N+1;
ENS = S*S';
ENEt = Et*Et';

norm = 1/ENEt;
Et = Et*norm;
Et = Et-mean(Et);

% nor = ENEt/ENS;
SN = Et;%S*nor;

z = zeros(1,2*N+1);

for k = 2:NN     % left
    p = k-1;
    z(k) = 0;
    for j = 1:p
        z(k) = z(k) + SN(p-j+1)*Et(N-j+1);
    end
end

for k = 2:N-1     % right
    z(N+k) = 0;
    for j = k:N
        z(N+k) = z(N+k) + SN(j)*Et(j-k+1);
    end
end

figure
subplot(211),plot(z),axis tight,grid on

for j = 1:length(z)
    u(j) = sum(z(1:j));
end

subplot(212),plot(u),axis tight,grid on