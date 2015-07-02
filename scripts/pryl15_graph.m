close all

%% pic 3
x = (1:Ts)/Fd/60;

figure('color','white')
subplot(211)
plot(x,val(1,:),'k'),grid,axis tight
ylabel('ЭКГ сигнал, мВ','FontName','Times New Roman','FontSize',12)
xlabel({'Время, мин' 'а) тестовый вариант обучающей выборки'},'FontName','Times New Roman','FontSize',12)

subplot(212)
plot(x,in_stab,'k'),grid,axis tight
ylabel('ЭКГ сигнал, мВ','FontName','Times New Roman','FontSize',12)
xlabel({'Время, мин' 'б) запись после стабилизация нулевой линии'},'FontName','Times New Roman','FontSize',12)

%% pic 4
x = (-47:80)/Fd*1e3;

figure('color','white')
plot(x,port{1},'k','LineWidth',1.5),hold on
plot(x,port{2},'--k','LineWidth',1.5),hold on
plot(x,port{3},':k','LineWidth',1.5),axis tight
legend('Normal beats','Atrial premature beats','Premature ventricular contraction')
xlabel('Время от R-зубца, мс'),ylabel('Форма QRS-комплексов'),grid

%% pic 5 Block 1
figure('color','white')
for i = 1:btypeN
   subplot(2,2,i+2),stem(des(i,:)*100,'.-k','LineWidth',1.5),axis([0 3 0 100]),grid
   set(gca,'XTickLabel','')
   xlabel(num2str(des(i,:)*100,3))
end
subplot(2,2,1:2),title({'Диаграммы качества распознавания QRS-комплексов по данным обучения ' ...
   ['(доля правильных решений – ' num2str((des(1,1)+des(2,2))/2e-2,2) '%)']},'FontName','Times New Roman','FontSize',12)
subplot(223),ylabel({'Частоты правильных и' 'ошибочных решений, %'},'FontName','Times New Roman','FontSize',12)

%% pic 6 Block 2
Btype = 1;
per = Bnum{Btype}(randi(length(Bnum{Btype}),1));
for btype = 1:btypeN
   cor_b2(btype) = nrm(f(per,:) - mean(f(per,:))) * nrm(port{btype} - mean(port{btype}))';
   cor_b2(btype) = (cor_b2(btype) +1)/2;
end
[~,Btype] = max(cor_b2);

figure('color','white')
stem(mean(cor,2)*100,'.-k','LineWidth',1.5),hold on
stem((1:2)+.1,cor_b2*100,'.:k','LineWidth',1.5),axis([0 3 0 100]),grid
legend('Среднее по данным обучения','Для текущего QRS-комплекса','FontName','Times New Roman','FontSize',12)
set(gca,'XTickLabel','')
ylabel({'Приведенные' 'корреляционные' 'интегралы'},'FontName','Times New Roman','FontSize',12)
title('Диаграмма качества распознавания QRS комплексов в текущей реализации','FontName','Times New Roman','FontSize',12)

%% pic 7
load('optwp')
[wid,pos] = find( optwp == max(max(optwp)) );
winwp = pos:pos+wid;
y = [min(port{3}) max(port{1})];
x = (-47:80)/Fd*1e3;
xwp = winwp/Fd*1e3+x(1)-1e3/Fd;

figure('color','white')
plot(x,port{1},'k',x,port{2},'--k',x,port{3},':k'),hold on
plot([xwp(end) xwp(end)],y,'--k','LineWidth',3),axis tight
plot(xwp,port{1}(winwp),'k','LineWidth',1.5),hold on
plot(xwp,port{2}(winwp),'--k','LineWidth',1.5),hold on
plot(xwp,port{3}(winwp),':k','LineWidth',1.5),hold on
plot([xwp(1) xwp(1)],y,'--k','LineWidth',3),hold on
legend('Normal beats','Atrial premature beats','Premature ventricular contraction','Borders of analysis window')
xlabel('????? ?? R-?????, ??'),ylabel({'?????????????? ?????' 'QRS-?????????? ?? ?????'}),grid

%% pic 8
load('indei')
x = (-47:80)/Fd*1e3;
num = 21;
xinden = setdiff(1:winL,wini{num});

figure('color','white')
plot(x,port{1},'-k',x,port{2},'--k',x,port{3},':k'),hold on
plot(xinden/Fd*1e3+x(1)-1e3/Fd,port{1}(xinden),'ok',xinden/Fd*1e3+x(1)-1e3/Fd,port{2}(xinden),'ok',...
   xinden/Fd*1e3+x(1)-1e3/Fd,port{3}(xinden),'ok'),axis tight
legend('Normal beats','Atrial premature beats','Premature ventricular contraction','??????????????? ??????')
xlabel('????? ?? R-?????, ??'),ylabel({'?????????????? ?????' 'QRS-?????????? ?? ?????'}),grid

%% pic 9
[~,maxi] = max(indei);
xinde = wini{maxi-1}/Fd*1e3+x(1)-1e3/Fd;

figure('color','white')
plot(x,port{1},'-k',x,port{2},'--k',x,port{3},':k'),hold on
plot(xinde,port{1}(wini{maxi-1}),'.k',xinde,port{2}(wini{maxi-1}),'.k',...
   xinde,port{3}(wini{maxi-1}),'.k'),axis tight
legend('Normal beats','Atrial premature beats','Premature ventricular contraction','????????????? ??????')
xlabel('????? ?? R-?????, ??'),grid
ylabel({'?????????????? ?????' 'QRS-?????????? ?? ?????'})

%% pic 10
figure('color','white')
plot(0:winL-1,indei*100,'.-k'),axis([0 winL-1 80 100]),grid
xlabel('?????????? ????????, ??????????? ?? ???????? ???????????????'),ylabel('???????? ?????????????, %')














