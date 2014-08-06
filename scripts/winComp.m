clc;
close all;
clear all;

load('Signals\Msasr')
n = 50;
s = 1;
input = [transform(Msasr(:,s:n:s+n*26-1),'vector') Msasr(1,s:n:s+n*26-1)];

N = floor(sqrt(length(input)));  % Chosen period
N1 = N^2;
Iter = 5;
PerLens = [26 20 13 7];
Per = length(PerLens);

Spec = zeros(Per,N);
for per = 1:Per
    disp(['period = ' num2str(PerLens(per))])
    
    data = input;
    Rep = zeros(Iter,PerLens(per));
    for iter = 1:Iter
        disp(['    iteration = ' num2str(iter)])

        W = win_OSR(data,PerLens(per));

        Rep(iter,:) = mean(W,1);

        lR = transform(Rep(iter,:),'vector_repeat');
        llR = transform(lR,'vector_repeat');
        dd = llR(1:length(data));
        data = data - dd;

    end
    
    lS = transform(sum(Rep,1),'vector_repeat');
    Spec(per,:) = lS(1:N);
end

% figure
% subplot(411),plot(Spec(1,:),'k'),axis tight,grid on,ylabel({'Длина';'периода';'26'})
% subplot(412),plot(Spec(2,:),'k'),axis tight,grid on,ylabel('20')
% subplot(413),plot(Spec(3,:),'k'),axis tight,grid on,ylabel('13')
% subplot(414),plot(Spec(4,:),'k'),axis tight,grid on,ylabel('7')

en = zeros(1,Per);
for i = 1:N-1
    en(i) = Spec(i,:)*Spec(i,:)';
end
figure
subplot(411),stem(en,'.k'),axis tight,grid on
    title('Распределение энергии','FontName','Times New Roman','FontSize',14)
    set(gca,'FontName','Times New Roman','FontSize',14,'XTick',2:26)

    
% figure
% subplot(511),plot(Rep(1,:),'k'),axis tight,grid on,ylabel({'№';'приближения';'1'})
% subplot(512),plot(sum(Rep(1:3,:),1),'k'),axis tight,grid on,ylabel('3')
% subplot(513),plot(sum(Rep(1:5,:),1),'k'),axis tight,grid on,ylabel('5')
% subplot(514),plot(sum(Rep(1:7,:),1),'k'),axis tight,grid on,ylabel('7')
% subplot(515),plot(sum(Rep(1:9,:),1),'k'),axis tight,grid on,ylabel('9')


% figure
% subplot(411),plot(W(1,:),'k'),axis tight,grid on
%     ylabel({'Номер';'окна';'';'';'1'},'FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     title('Периодические составляющие','FontName','Times New Roman','FontSize',14)
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% subplot(412),plot(W(11,:),'k'),axis tight,grid on
%     ylabel('11','FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     title('. . . . . . . . . . . . . . . . . . . . . . . . . . .','FontSize',24,'FontWeight','Bold')
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% subplot(413),plot(W(12,:),'k'),axis tight,grid on
%     ylabel('12','FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     title('. . . . . . . . . . . . . . . . . . . . . . . . . . .','FontSize',24,'FontWeight','Bold')
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% subplot(414),plot(W(13,:),'k'),axis tight,grid on
%     ylabel('13','FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% figure
% subplot(411),plot(W(14,:),'k'),axis tight,grid on
%     ylabel('14','FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% subplot(412),plot(W(17,:),'k'),axis tight,grid on
%     ylabel('17','FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     title('. . . . . . . . . . . . . . . . . . . . . . . . . . .','FontSize',24,'FontWeight','Bold')
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% subplot(413),plot(W(27,:),'k'),axis tight,grid on
%     ylabel('27','FontName','Times New Roman','FontSize',14,'Rotation',0,'FontWeight','Bold')
%     title('. . . . . . . . . . . . . . . . . . . . . . . . . . .','FontSize',24,'FontWeight','Bold')
%     set(gca,'FontName','Times New Roman','FontSize',14,'XTick',1:26)
% 
% figure,surf(W)