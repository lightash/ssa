clc;
close all;
clear all;

load('Signals\Msasr')
N = 26;
s = 1;%503;
n = 50;%1;
a1 = Msasr(1:N,s:n:s+n*N-1);
N1 = N*N;

srednee = mean(a1);
ss = transform(srednee,'vector_repeat');

% correlator on N periods
Et = ss;
shum = .1*randn(1,N1);
sdvig = 0;
S = [zeros(1,sdvig) ss(1:length(ss)-sdvig)] + shum;

z = zeros(N-1 + 1,N);
for i = 1:N-1
    for j = 1:N
        k = (i-1)*N+1 + j-1 : i*N + j-1;
        ENEt = sqrt(Et(k)*Et(k)');
        ENS = sqrt(S(k)*S(k)');
        z(i,j) = ( S(k)/ENS )*( Et(k)/ENEt )';
    end
    z(N,:) = z(N,:) + z(i,:);
end

figure
subplot(311),plot(S),axis tight,grid on
subplot(312),plot(z(1:N-1,:)'),axis tight,grid on
subplot(313),plot(z(N,:)),axis tight,grid on

disp(sum(z(N,:)))
