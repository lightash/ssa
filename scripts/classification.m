clc;
close all;
clear all;

load('Signals\krov','krov2','krov3')
N = 21;
N1 = N^2;
s = 148;

k11 = krov2( s : s+N1-1, 2);
k12 = krov2( s+N1 : s+2*N1-1, 2);
k21 = krov3( s : s+N1-1, 2);
k22 = krov3( s+N1 : s+2*N1-1, 2);

figure('Color','w')
subplot(211),plot(krov2(:,2),'k'),axis tight,grid on
    title('»сходна€ реализаци€ дл€ объекта класса 1','FontName','Times New Roman','FontSize',14)
subplot(212),plot(krov3(:,2),'k'),axis tight,grid on
    title('»сходна€ реализаци€ дл€ объекта класса 2','FontName','Times New Roman','FontSize',14)

A1_11 = transform(k11,'matrix');
A1_12 = transform(k12,'matrix');
A1_21 = transform(k21,'matrix');
A1_22 = transform(k22,'matrix');

B_11 = imp_OSR(A1_11);
B_12 = imp_OSR(A1_12);
B_21 = imp_OSR(A1_21);
B_22 = imp_OSR(A1_22);

Bs_11 = sum(B_11,2);
Bs_12 = sum(B_12,2);
Bs_21 = sum(B_21,2);
Bs_22 = sum(B_22,2);

Ebs_11 = zeros(1,21);
Ebs_12 = zeros(1,21);
Ebs_21 = zeros(1,21);
Ebs_22 = zeros(1,21);
for i = 1:21
    Ebs_11(i) = B_11(:,i)'*B_11(:,i);
    Ebs_12(i) = B_12(:,i)'*B_12(:,i);
    Ebs_21(i) = B_21(:,i)'*B_21(:,i);
    Ebs_22(i) = B_22(:,i)'*B_22(:,i);
end

% ќбъект класса 1   ќбъект класса 2
% ѕериодическа€ составл€юща€ фрагмента 1
% Ёнергетический спектр фрагмента 1
figure('Color','w') 
subplot(421),plot(Bs_11,'k'),axis tight,grid on
    title({'ќбъект класса 1';'ѕериодическа€ составл€юща€ фрагмента 1'},'FontName','Times New Roman','FontSize',14)
subplot(422),plot(Bs_21,'k'),axis tight,grid on
    title({'ќбъект класса 2';'ѕериодическа€ составл€юща€ фрагмента 1'},'FontName','Times New Roman','FontSize',14)
subplot(423),stem(Ebs_11,'.k'),axis tight,grid on
    title('Ёнергетический спектр фрагмента 1','FontName','Times New Roman','FontSize',14)
subplot(424),stem(Ebs_21,'.k'),axis tight,grid on
    title('Ёнергетический спектр фрагмента 1','FontName','Times New Roman','FontSize',14)
subplot(425),plot(Bs_12,'k'),axis tight,grid on
    title('ѕериодическа€ составл€юща€ фрагмента 2','FontName','Times New Roman','FontSize',14)
subplot(426),plot(Bs_22,'k'),axis tight,grid on
    title('ѕериодическа€ составл€юща€ фрагмента 2','FontName','Times New Roman','FontSize',14)
subplot(427),stem(Ebs_12,'.k'),axis tight,grid on
    title('Ёнергетический спектр фрагмента 2','FontName','Times New Roman','FontSize',14)
subplot(428),stem(Ebs_22,'.k'),axis tight,grid on
    title('Ёнергетический спектр фрагмента 2','FontName','Times New Roman','FontSize',14)

% figure('Color','w')
% subplot(221),plot(k11,'k'),axis tight,grid on
%     title({'–еализаци€ дл€ объекта класса 1,';'фрагмент 1'},'FontName','Times New Roman','FontSize',14)
% subplot(222),plot(k12,'k'),axis tight,grid on
%     title({'–еализаци€ дл€ объекта класса 1,';'фрагмент 2'},'FontName','Times New Roman','FontSize',14)
% subplot(223),plot(k21,'k'),axis tight,grid on
%     title({'–еализаци€ дл€ объекта класса 2,';'фрагмент 1'},'FontName','Times New Roman','FontSize',14)
% subplot(224),plot(k22,'k'),axis tight,grid on
%     title({'–еализаци€ дл€ объекта класса 1,';'фрагмент 2'},'FontName','Times New Roman','FontSize',14)

% data = B_11;
% 
% axes = zeros(N,N);
% proj = axes;
% vproj = zeros(N);
% for axi = 1:N
%     
%     axes(axi,:) = data(axi,:);
% %     axis(axi,:) = axis(axi,:)/sqrt(axis(axi,:) * axis(axi,:)');
%     
%     for i = 1:N
%         proj(axi,i) = data(i,:) * axes(axi,:)';  % Lenghts of projetions on chosen axes
%         vproj(i,:) = proj(axi,i) * axes(axi,:);  % Vectors of projections
%     end
%     
%     data = data - vproj;
%     
% end
% 
% for axi = 1:N
%     vproj(axi,:) = proj(axi,i)*axes(axi,:);
% end
% 
% figure
% subplot(421),plot(B_11(1,:)),axis tight,grid on
% subplot(422),plot(vproj(1,:)),axis tight,grid on
% subplot(423),plot(B_11(2,:)),axis tight,grid on
% subplot(424),plot(vproj(2,:)),axis tight,grid on
% subplot(425),plot(B_11(20,:)),axis tight,grid on
% subplot(426),plot(vproj(20,:)),axis tight,grid on
% subplot(427),plot(B_11(21,:)),axis tight,grid on
% subplot(428),plot(vproj(21,:)),axis tight,grid on
% 
% 
% Ebs_11 = zeros(1,21);
% for i = 1:21
%     Ebs_11(i) = vproj(i,:)*vproj(i,:)';
% end
% 
% figure
% plot(Ebs_11),axis tight,grid on
