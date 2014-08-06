clc;
close all;
clear all;

load('Signals\Msasr')
N = 26;
s = 1;%503;
n = 50;%1;
a1 = Msasr(1:N,s:n:s+n*N-1);
N1 = N*N;

ss = transform(a1,'vector');

B = imp_OSR(a1);
Et = transform(sum(B,2)','vector_repeat');

% filter with 1 period
s = Et;                 % эталон

shum = 1*randn(1,N1);
u = ss + shum;          % вход
% % N = N1;

umat = transform(u,'matrix');
smat = transform(s,'matrix');
for i = 1:N
    ENu = umat(i,:)*umat(i,:)';
    umat(i,:) = umat(i,:)/sqrt(ENu);
    ENs = smat(i,:)*smat(i,:)';
    smat(i,:) = smat(i,:)/sqrt(ENs);
end
u = transform(umat,'vector');
s = transform(smat,'vector');
% u = u*sqrt(26);

h = zeros(1,N1);
for i = 1:N1
    h(i) = s(N1-i+1);  % имп. хар-ка
end
% h = h*sqrt(26);

t0 = N1;

uv = zeros(1,N1);
sv = zeros(1,N1);
st = zeros(1,N1);
for t = 1:N1
    for x = 1:t
        uv(t) = uv(t) + u(x)*h(t-x+1);  % выход
    end
end

figure('Color','w')
subplot(311),plot(u,'k'),axis tight,grid on
    title('Исходный сигнал','FontName','Times New Roman','FontSize',14)
subplot(312),plot(h(1:N),'k'),axis tight,grid on
    title('Импульсная характеристика СФ','FontName','Times New Roman','FontSize',14)
subplot(313),plot(uv,'k'),axis tight,grid on
    title('Выходной сигнал СФ','FontName','Times New Roman','FontSize',14)

