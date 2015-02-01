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

% for k = 1:1
    shum = .0*randn(1,N1);
    S = ss + shum;

    z = zeros(N+1,2*N-1);
    for i = 1:N
        sig(i,:) = S( (i-1)*N+1 : i*N );
        ENsig = sig(i,:)*sig(i,:)';
        sig(i,:) = sig(i,:)/sqrt(ENsig);
        
        etal(i,:) = Et( (i-1)*N+1 : i*N );
        ENetal = etal(i,:)*etal(i,:)';
        etal(i,:) = etal(i,:)/sqrt(ENetal);
        
        z(i,:) = xcorr(etal(i,:),sig(i,:));
        z(N+1,:) = z(N+1,:) + z(i,:);
    end

%     fmax(k) = z(N+1,26);
%     smax(k) = max([z(N+1,1:25) z(N+1,27:51)]);
% end

figure('Color','w')
subplot(411),plot(S,'k'),axis tight,grid on
    title('Исходный сигнал','FontName','Times New Roman','FontSize',14)
subplot(412),plot(Et(1:N),'k'),axis tight,grid on
    title('Эталон для одного периода','FontName','Times New Roman','FontSize',14)
subplot(413),plot(z(1:N,:)','k'),axis tight,grid on
    title('Корреляционные функции по периодам','FontName','Times New Roman','FontSize',14)
subplot(414),plot(z(N+1,:),'k'),axis tight,grid on
    title('Суммарная корреляционная функция периодов','FontName','Times New Roman','FontSize',14)

% figure
% plot(fmax,'b'),hold on,plot(smax,'r'),axis tight,grid on
