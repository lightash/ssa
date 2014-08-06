clc;
close all;
clear all;

% For section 2
load('Signals\Msasr')
n = 50;
s = 1;
input = transform(Msasr(:,s:n:s+n*26-1),'vector');

N = floor(sqrt(length(input)));  % Chosen period
N1 = N^2;

data = transform(input,'matrix');

B = imp_OSR(data);

figure
subplot(311),plot(B(:,2),'k'),axis tight,grid on
    title('»мпульсные составл€ющие'),ylabel('1')
subplot(312),plot(B(:,3),'k'),axis tight,grid on,ylabel('2')
subplot(313),plot(B(:,4),'k'),axis tight,grid on,ylabel('3')

disp(B(:,2)'*B(:,2))
disp(B(:,3)'*B(:,3))
disp(B(:,4)'*B(:,4))
disp(' ')

figure
subplot(311),plot(B(:,24),'k'),axis tight,grid on,ylabel('23')
subplot(312),plot(B(:,25),'k'),axis tight,grid on,ylabel('24')
subplot(313),plot(B(:,26),'k'),axis tight,grid on,ylabel('25')

disp(B(:,24)'*B(:,24))
disp(B(:,25)'*B(:,25))
disp(B(:,26)'*B(:,26))
disp(' ')

en = zeros(1,N-1);
for i = 2:N
    en(i-1) = B(:,i)'*B(:,i);
end

figure
subplot(311),plot(sum(B(:,2:end),2),'k'),axis tight,grid on
    title('—умма импульсных составл€ющих')
subplot(312),stem(en,'.k'),axis tight,grid on
    title('–аспределение энергии')

disp(sum(B(:,2:end),2)'*sum(B(:,2:end),2))
disp(' ')


figure
subplot(311),plot(sum(B,2),'k'),axis tight,grid on
    title('—умма нерперывной и импульсных составл€ющих')
subplot(312),plot(input(1:N),'k'),axis tight,grid on
    title('»сходный сигнал')

disp(sum(B,2)'*sum(B,2))
disp(input(1:N)*input(1:N)')
disp(' ')



% For section 3
load('Signals\krov','krov2')
N = 21;
N1 = N^2;
s = 148;

k = krov2( s : s+N1-1, 2);
A1 = transform(k,'matrix');
[B, A3] = imp_OSR(A1);

Bs = sum(B,2);
Ebs = zeros(1,21);
for i = 1:21
    Ebs(i) = B(:,i)'*B(:,i);
end

figure('Color','w')
subplot(411),plot(k,'k'),axis tight,grid on
    title('»сходный сигнал','FontName','Times New Roman','FontSize',14)
subplot(412),plot(B(:,1),'k'),axis tight,grid on
    title('"Ќепрерывна€" составл€юща€','FontName','Times New Roman','FontSize',14)
subplot(413),plot(B(:,2),'k'),axis tight,grid on
    title('»мпульсные составл€ющие','FontName','Times New Roman','FontSize',14)
    ylabel('1','FontName','Times New Roman','FontSize',14,'Rotation',0)
subplot(414),plot(B(:,3),'k'),axis tight,grid on
    ylabel('2','FontName','Times New Roman','FontSize',14,'Rotation',0)
figure('Color','w')
str = '. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .';
subplot(411),plot(B(:,20),'k'),axis tight,grid on
    title(str,'FontName','Times New Roman','FontSize',24)
    ylabel('19','FontName','Times New Roman','FontSize',14,'Rotation',0)
subplot(412),plot(B(:,21),'k'),axis tight,grid on
    ylabel('20','FontName','Times New Roman','FontSize',14,'Rotation',0)
subplot(413),stem(Ebs,'.k'),axis tight,grid on
    title('Ёнергетический спектр составл€ющих','FontName','Times New Roman','FontSize',14)
figure('Color','w')
subplot(411),plot(Bs,'k'),axis tight,grid on
    title('—уммарна€ периодическа€ составл€юща€','FontName','Times New Roman','FontSize',14)
subplot(412),plot(transform(A3,'vector'),'k'),axis tight,grid on
    title('ќстатки аппроксимации','FontName','Times New Roman','FontSize',14)

    