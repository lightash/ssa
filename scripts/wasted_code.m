tt0 = 2*pi;
tt1 = tt0/N1;
tt2 = tt1/16;
tt3 = tt2/4;
te1 = tet1;
te2 = tet2/16;
te3 = tet3/64;
dp = mapi/16;

j2 = 0;
j3 = 0;
j4 = 0;
as1 = zeros(1,N1);
as2 = zeros(1,N1);
api = zeros(1,N1);
as3 = zeros(1,N1);
ame = zeros(1,N1);
a3 = zeros(1,N1);
S = zeros(1,N1);

for i = 1:N1
   as1(i) = mas1*sin( te1 + (i-1)*tt1 );
   j2 = j2+1;
   if j2 == 17
      j2 = 1;
      api(i) = -mapi/2 + (j2-1)*dp;         % mapi2/2 - ?
   end
   as2(i) = mas2*sin( te2 + (j2-1)*tt2 );
   j3 = j3+1;
   if j3 == 5
       j3 = 1;
   end
   as3(i) = mas3*sin( te3 + (j3-1)*tt3 );
   if j3 == 1 || j3 == 2
       ame(i) = mame;
   elseif j3 == 3 || j3 == 4
       ame(i) = -mame;
   end
   j4 = j4+1;
   if j4 == 3
      j4 = 1;
   end
   if j4 == 1
       a3(i) = a31;
   elseif j4 == 2
       a3(i) = a32;
   end
   S(i) = as1(i) + as2(i) + api(i) + as3(i) + ame(i) + a3(i) + apost;
end

ms = 0;
for i = 1:N1
    ms = ms + S(i)*S(i);
end
ms = sqrt(ms);
ee = S./ms;

%------------------------------------#9------------------------------------
% Формирование матрицы DEL[1:N1-1][1:N1-1]

DEL = zeros(N1-1);

for i = 1:N1-1
   for j = 1:N1-1
       DEL(i,j) = a2(j,i);
   end
end

%------------------------------------#10-----------------------------------
% Формирование вектора G[1:N1-1]

G = zeros(1,N1-1);

for j = 1:N1-1
   G(j) = a2(N1,j); 
end

%------------------------------------#11-----------------------------------
% Расчет K[1:N1-1] из СЛАУ ||DEL|| * ||K|| = ||G||

K = DEL\(G');

sha1 = zeros(N1);
i = 1;
j = 1;
for k = 1:N1*N1
    sha1(i,j) = sh(k);
    if j < N1
        j = j+1;
    else
        j = 1;
        i = i+1;
    end
end

% --- матрица -> вектор ---
ssa = sha2(1,:);
for i = 2:N1
    ssa = [ssa sha2(i,:)];
end
sig = S + ssa;

% --- вектор -> матрица ---
a1 = zeros(N1);
i = 1;
j = 1;
for k = 1:N1*N1
    a1(i,j) = sig(k);
    if j < N1
        j = j+1;
    else
        j = 1;
        i = i+1;
    end
end

% --- вектор -> матрица ---
sha1(1,:) = sh(1:N1);
for i = 1:N1
   sha1(i,:) = sh( (i-1)*N1+1 : i*N1 );
end

%--------------------------------------------------------------------------
[sha2, shb2] = OSR(sha1, N1);
%--------------------------------------------------------------------------

shSK = shb2*sha2';
disp(['lg(sh) -> ' num2str(round(log10(norm(shSK,inf))))])

--- вектор -> матрица ---
sh1(1,:) = sh(1:N1);
for i = 1:N1
   sh1(i,:) = sh( (i-1)*N1+1 : i*N1 );
end

%------------------------------------#12-----------------------------------
% Расчет (b2,a2) - скалярного произведения для проверки ортогональности

SK = b2*a2';s
disp(['lg(SK)_' num2str(N1) ' -> ' num2str(round(log10(norm(SK,inf))))])

t = 0:1/N1:1-1/N1;
si = sin(2*pi*t);
a1 = a1 + [si ; si ; si ; si];

N1 = 8;
a1 = [a1 a1 ; a1 a1];

% Вычитка в .mat из .txt изменения дифференциальной температуры воды с
%   кровью разных концентраций
clc;close all;clear all;
file1 = fopen('D:\Политех\Диплом\krov\26.09.2012 11\krov1.txt');
file2 = fopen('D:\Политех\Диплом\krov\26.09.2012 22\krov-2.txt');
file3 = fopen('D:\Политех\Диплом\krov\26.09.2012 33\krov 3.txt');
br=0;
for k = 1:inf
    str1 = fgetl(file1);
    if ischar(str1)==1
        krov1(k,:) = sscanf(str1,'%f%*s%f',[1,inf]);
    else
        br=br+1;
    end
    str2 = fgetl(file2);
    if ischar(str2)==1
        krov2(k,:) = sscanf(str2,'%f%*s%f',[1,inf]);
    else
        br=br+1;
    end
    str3 = fgetl(file3);
    if ischar(str3)==1
        krov3(k,:) = sscanf(str3,'%f%*s%f',[1,inf]);
    else
        br=br+1;
    end
    if br==3
        break
    end
end
save('D:\Dropbox\MATLAB\Диплом\krov.mat','krov1','krov2','krov3')

% Что-то про корреляцию
xc(1,:) = xcorr(a1(1,:),b2);%[b2 zeros(1,12)]);
xc(2,:) = xcorr(a1(2,:),b2);%[zeros(1,4) b2 zeros(1,8)]);
xc(3,:) = xcorr(a1(3,:),b2);%[zeros(1,8) b2 zeros(1,4)]);
xc(4,:) = xcorr(a1(4,:),b2);%[zeros(1,12) b2]);
xcs = xc(1,:)+xc(2,:)+xc(3,:)+xc(4,:);%xc(1,12:31)+xc(2,8:27)+xc(3,4:23)+xc(4,1:20);
co = 'k';
figure%hold on
subplot(511),plot(xc(1,:),co),axis tight,ylabel('a1'),grid
subplot(512),plot(xc(2,:),co),axis tight,ylabel('a2'),grid
subplot(513),plot(xc(3,:),co),axis tight,ylabel('b2'),grid
subplot(514),plot(xc(4,:),co),axis tight,ylabel('correlation'),grid
subplot(515),plot(xcs,co),axis tight

a11 = [b2 ; b2 ; b2 ; b2];
xc(1,:) = xcorr(a11(1,:),b2);
xc(2,:) = xcorr(a11(2,:),b2);
xc(3,:) = xcorr(a11(3,:),b2);
xc(4,:) = xcorr(a11(4,:),b2);
xcs = xc(1,:)+xc(2,:)+xc(3,:)+xc(4,:);

a11n = a11 + 4*randn(4);
xcn(1,:) = xcorr(a11n(1,:),b2);
xcn(2,:) = xcorr(a11n(2,:),b2);
xcn(3,:) = xcorr(a11n(3,:),b2);
xcn(4,:) = xcorr(a11n(4,:),b2);
xcsn = xcn(1,:)+xcn(2,:)+xcn(3,:)+xcn(4,:);

co = '--';
figure
subplot(511),plot(1:2*N1-1,xc(1,:),1:2*N1-1,xcn(1,:),co),grid,axis tight
subplot(512),plot(1:2*N1-1,xc(2,:),1:2*N1-1,xcn(2,:),co),grid,axis tight
subplot(513),plot(1:2*N1-1,xc(3,:),1:2*N1-1,xcn(3,:),co),grid,axis tight
subplot(514),plot(1:2*N1-1,xc(4,:),1:2*N1-1,xcn(4,:),co),grid,axis tight
subplot(515),plot(1:2*N1-1,xcs,1:2*N1-1,xcsn,co),grid,axis tight

% Какие-то итерации
bb(1,1:N1) = b2;
S(1,1:N1) = mean(a2);
a2ost = zeros(N1);
for i = 1:N1
    a2ost(i,:) = a2(i,:) - S(1,1:N1);
end
a2st(1,:) = a2ost(N1,:);
a2st(2,:) = a2ost(:,N1)';

N1 = N1-1;
a2min = a2ost(1:N1,1:N1);

% 2

a2min = preprocess(a2min, 0);
[a3, b3] = OSR(a2min);

bb(2,1:N1) = b3;
S(2,1:N1) = mean(a3);
a3ost = zeros(N1);
for i = 1:N1
    a3ost(i,:) = a3(i,:) - S(2,1:N1);
end
a3st(1,:) = a3ost(N1,:);
a3st(2,:) = a3ost(:,N1)';

N1 = N1-1;
a3min = a3ost(1:N1,1:N1);

% 3

a3min = preprocess(a3min, 0);
[a4, b4] = OSR(a3min);

bb(3,1:N1) = b4;
S(3,1:N1) = mean(a4);
a4ost = zeros(N1);
for i = 1:N1
    a4ost(i,:) = a4(i,:) - S(3,1:N1);
end
a4st(1,:) = a4ost(N1,:);
a4st(2,:) = a4ost(:,N1)';


figure
subplot(311),plot(a1(1,:)),axis tight,ylabel('a1_1'),grid
subplot(312),plot(a2min(1,:)),axis tight,ylabel('a2min_1'),grid
subplot(313),plot(a3min(1,:)),axis tight,ylabel('a3min_1'),grid

figure
subplot(311),plot(a2(1,:)),axis tight,ylabel('a2_1'),grid
subplot(312),plot(a3(1,:)),axis tight,ylabel('a3_1'),grid
subplot(313),plot(a4(1,:)),axis tight,ylabel('a4_1'),grid

figure
subplot(311),plot(bb(1,:)),axis tight,ylabel('bb_1'),grid
subplot(312),plot(bb(2,:)),axis tight,ylabel('bb_2'),grid
subplot(313),plot(bb(3,:)),axis tight,ylabel('bb_3'),grid

figure
subplot(311),plot(S(1,:)),axis tight,ylabel('S_1'),grid
subplot(312),plot(S(2,:)),axis tight,ylabel('S_2'),grid
subplot(313),plot(S(3,:)),axis tight,ylabel('S_3'),grid

figure
subplot(311),plot(a2ost(1,:)),axis tight,ylabel('a2ost_1'),grid
subplot(312),plot(a3ost(1,:)),axis tight,ylabel('a3ost_1'),grid
subplot(313),plot(a4ost(1,:)),axis tight,ylabel('a4ost_1'),grid

% Увеличение размера матрицы в 2 раза
N1 = sqrt(N1);
a12 = b2(1:N1);
for i = 2:N1
    a12 = [a12 ; b2((i-1)*N1+1:i*N1)];
end

a12 = preprocess(a12, 0);
[a22, b22] = OSR(a12);

a = a22(1,:); b = b22; aa = a12(1,:);
for i = 2:N1
    a = [a a22(i,:)];
    b = [b  b22];
    aa = [aa a12(i,:)];
end

figure
subplot(311),plot(aa),axis tight,ylabel('a1'),grid
subplot(312),plot(a),axis tight,ylabel('a2'),grid
subplot(313),plot(b),axis tight,ylabel('b2'),grid


N1 = sqrt(N1);
a13 = [b22(1:2) ; b22(3:4)];

a13 = preprocess(a13, 0);
[a23, b23] = OSR(a13);

a = a23(1,:); b = [b23 b23]; aa = a13(1,:);
for i = 2:N1
    a = [a a23(i,:)];
    b = [b  b23 b23];
    aa = [aa a13(i,:)];
end

figure
subplot(311),plot(aa),axis tight,ylabel('a1'),grid
subplot(312),plot(a),axis tight,ylabel('a2'),grid
subplot(313),plot(b),axis tight,ylabel('b2'),grid


b24 = [ abs( b23(2)-b23(1) )/2  -abs( b23(2)-b23(1) )/2 ];
b25 = ( b23(1)+b23(2) )/2;

b = b24; bb = [b25 b25];
for i = 2:8
    b = [b b24];
    bb = [bb b25 b25];
end

figure
subplot(211),plot(b),grid,axis tight,ylabel('b24')
subplot(212),plot(bb),grid,axis tight,ylabel('b25')


b22_och = b22-[b23 b23];
b2_och = b2-[b22_och b22_och b22_och b22_och]-[b23 b23 b23 b23 b23 b23 b23 b23];

figure
subplot(211),plot([b22_och b22_och b22_och b22_och]),grid,axis tight,ylabel('b22_och')
subplot(212),plot(b2_och),grid,axis tight,ylabel('b2_och')

w1 = 0;
for i = 1:16
    w1 = w1 + b25*b25;
end

w2 = [b24 b24 b24 b24 b24 b24 b24 b24]*[b24 b24 b24 b24 b24 b24 b24 b24]';

w4 = [b22_och b22_och]*[b22_och b22_och]';

w16 = b2_och*b2_och';

W(1) = w1;
W(2) = w2;
W(3) = w4;
W(4) = w16;
figure,stem(log10(W)),grid,axis tight,title(W)

% Уменьшение размера матрицы на 1
N1 = N1-1;

mod_b2n = zeros(N1+1);

for z1 = 1:N1+1
    for z2 = 1:N1+1
        a1n = [a1(1:z1-1,1:z2-1) a1(1:z1-1,z2+1:N1+1)];
        a1n = [a1n ; a1(z1+1:N1+1,1:z2-1) a1(z1+1:N1+1,z2+1:N1+1)];
        [a2n, b2n] = OSR(a1n, N1);
        mod_b2n(z1,z2) = sqrt(b2n*b2n');
    end
end
max_mod_b2n = max(max(mod_b2n));
[i_z1 i_z2] = find(mod_b2n == max_mod_b2n);
disp(['z1=' num2str(i_z1) '; z2=' num2str(i_z2) '; max=' num2str(max_mod_b2n)])
figure,surf(mod_b2n),axis tight,xlabel('z2'),ylabel('z1')
title(['z1=' num2str(i_z1) '; z2=' num2str(i_z2) '; max=' num2str(max_mod_b2n)])

%% Оконная обработка
sr = mean(a1);           % Среднее значение строк
mod_e = sqrt(sr*sr');  	% Длина орта среднего
e_sr = sr./mod_e;        % Орт среднего

a = a1(1,:);
for i = 2:N1
    a = [a a1(i,:)];    % Разложение a0 в строку a
end

SSCh = zeros(1,2*N1*N1+1);

for k = 1:20                        % По наборам случайных чисел
    disp([num2str(k/20*100) '%'])
    
    % Добавление шума по бокам
    sig = [5*randn(1,N1*N1) a+5*randn(1,N1*N1) 5*randn(1,N1*N1)];

    SSC = zeros(1,2*N1*N1+1);
    for z = 1:2*N1*N1+1             % Проход окном по сигналу
        okno = sig(z:z-1+N1*N1);

        % Окно -> матрица
        a1 = zeros(N1);
        for i = 1:N1
            a1(i,:) = okno((i-1)*N1+1:(i-1)*N1+N1);     % Матрица окна
        end

        % Проецирование исходного сигнала на орт среднего
        a_nor = zeros(N1);
        for i = 1:N1
            a_nor(i,:) = a0(i,:)*(e_sr') * e_sr;
        end

        a_sr = mean(a_nor);
        ma_sr = a_sr*a_sr';     % ~sqrt?

        p_nor = zeros(1,N1);
        for i = 1:N1
            p_nor(i) = sqrt(a_nor(i,:)*(a_nor(i,:)')) / ma_sr;  % ?
        end

        a1_nor = zeros(N1);
        for i = 1:N1
            a1_nor(i,:) = a1(i,:) / p_nor(i);
        end

        a1_dop_nor = a1_nor / ma_sr;

        SC = zeros(1,N1);
        for i = 1:N1
            SC(i) = a1_dop_nor(i,:)*e_sr';
        end

        SSC(z) = sum(SC);
    end
    SSCh = SSCh + SSC;
end
SSCh = SSCh/10;
figure,plot(SSCh)

%% Графики из книжки
N1 = 4;
blood1 = [1.47 1.08 .77 .53  .37 .32 .31 .33  .35 .37 .4  .35  .3 .26 .2 .16];
blood2 = [1.2 1.03 .88 .78  .72 .71 .8 .93  1.03 1  .9 .76  .68 .47 .3 .21];
figure,plot(1:N1^2,blood1,1:N1^2,blood2,'r'),grid,axis tight
a1 = blood1;

%% Тестовый сигнал
N1 = 256;
mas1 = 1;       % Амплитуды синусоид
mas2 = 0.5;
mas3 = 0.25;
mapi = 0.75;    % Амплитуда пилы
mame = 0.35;    %           меандра
a31 = 0.1;      % Повторяющаяся пара чисел
a32 = 0.3;
apost = 1;      % Постоянная составляющая
tet1 = 0;       % Начальные фазы синусоид
tet2 = 0;
tet3 = 0;

S = testSignal(N1,[mas1,mas2,mas3,mapi,mame,a31,a32,apost],[tet1,tet2,tet3]);

%% Новый базис
bas = [1 1 ; -1 1];
bas = [bas bas ; bas -1*bas];
bas = [bas bas ; bas -1*bas];
bas = [bas bas ; bas -1*bas];

for i = 1:N1
    a1(i) = a1(i,:)*bas(1,:)'/sqrt(N1);
end

mean = mean(a1);

sr = mean*bas(1,:)./sqrt(N1);

for i = 1:N1
    a11(i,:) = a1(i,:)-sr;
end


for j = 1:N1
    for i = 1:N1
        if i <= 6 || i == 15 || i == 16
            bask(j,i,:) = bas(i,:)*a11(j,i);
        else
            bask(j,i,:) = 0;
        end
    end
end
a111 = zeros(N1);
for j = 1:N1
    for i = 1:N1
        string(1,:) = bask(j,i,:);
        a111(j,:) = a111(j,:) + string;
    end
end
a111_ost = a1 - a111;

figure,plot(a1')
figure,plot(a111')
figure,plot(a111_ost')

%% Внутренние составляющие
a1_1 = [b2(1:2);b2(3:4)];
[a2_1,b2_2] = OSR(a1_1);
b2_3 = mean(b2_2);

b_1 = b2_3;
for i = 2:16
    b_1 = [b_1 b2_3];
end
b_2 = b2_2;
for i = 2:8
    b_2 = [b_2 b2_2];
end
b_4 = [b2 b2 b2 b2];

b1_4 = b_4 - b_2 - b_1;
b1_2 = b_2 - b_1;

figure
subplot(311),plot(b_4,'.--k'),axis tight,ylabel('b_4'),grid
subplot(312),plot(b_2,'.--k'),axis tight,ylabel('b_2'),grid
subplot(313),plot(b_1,'.--k'),axis tight,ylabel('b_1'),grid
% subplot(514),plot(b1_4,'.--k'),axis tight,ylabel('b''_4'),grid,legend('b''_4=b_4-b_2-b_1')
% subplot(515),plot(b1_2,'.--k'),axis tight,ylabel('b''_2'),grid,legend('b''_2=b_2-b_1')

%%
a1_1 = a2(1:3,1:3);
[a2_1,b2_1] = OSR(a1_1);
b2_1 = [b2_1 0];
b2_1 = [b2_1 b2_1 b2_1 b2_1];
b2s(1,:) = b(1:N1)+b2_1(1:N1);
for i = 1:N1
   b2s(i,:) = b( (i-1)*N1+1 : i*N1 )+b2_1( (i-1)*N1+1 : i*N1 );
end
[a2_2,b2_2] = OSR(a1-b2s);
a2s = a2_2(1,:);
for i = 2:N1
    a2s = [a2s a2_2(i,:)];
end

figure
subplot(411),plot(b2_1),axis tight,grid on
subplot(412),plot(a-b2_1),axis tight,grid on
subplot(413),plot(b+b2_1),axis tight,grid on
subplot(414),plot(a2s),axis tight,grid on

a1_2 = a2_2(1:2,1:2);
[a2_3,b2_3] = OSR(a1_2);
b2_3 = [b2_3 0 0];
b2_3 = [b2_3 b2_3 b2_3 b2_3];

figure,plot(b2_3+b+b2_1),axis tight,grid on

b3s(1,:) = b(1:N1)+b2_1(1:N1)+b2_3(1:N1);
for i = 1:N1
   b3s(i,:) = b( (i-1)*N1+1 : i*N1 )+b2_1( (i-1)*N1+1 : i*N1 )+b2_3( (i-1)*N1+1 : i*N1 );
end
[a2_4,b2_4] = OSR(a1-b3s);
a3s = a2_4(1,:);
for i = 2:N1
    a3s = [a3s a2_4(i,:)];
end
figure,plot(a3s),axis tight, grid on

%% innerComp
load('Signals\Msasr')
N1 = 26;
s = 1;%503;
n = 50;%1;
a1 = Msasr(1:N1,s:n:s+n*N1-1);

% [a1] = preprocess(a1, 0);
[a2, b2, e2, p] = OSR(a1);

% Построить всё
a = a2(1,:);b = b2;aa = a1(1,:);
for i = 2:N1
    a = [a a2(i,:)];
    b = [b  b2];
    aa = [aa a1(i,:)];
end
% figure
% subplot(311),plot(aa),axis tight,title('Исходная последовательность a1'),grid
% subplot(312),plot(b),axis tight,title('Периодическая составляющая b2'),grid
% subplot(313),plot(a),axis tight,title('Остатки a2'),grid

A1 = zeros(N1,N1,N1);
A1(:,:,N1) = a1;
A2 = zeros(N1,N1,N1);
A2(:,:,N1) = a2;
B2 = zeros(1,N1,N1);
B2(1,:,N1) = b2;

IJ = 1:N1;
figure
for size = N1-1:-1:N1-4
    disp(size)
    
    Eb2 = zeros(1,size);
    ijr = 0;
    for k = 1:size
        ijr = IJ(1:k-1);
        ijr = [ijr IJ(k+1:length(IJ))];
        ar = A2(ijr,ijr,size+1);
%         [ar] = preprocess(ar, 0);
        [~, b2] = OSR(ar);
        Eb2(k) = b2*b2';
    end
    [~,ijm] = max(Eb2);
    i = IJ(1:ijm-1);
    IJ = [i IJ(ijm+1:length(IJ))];
    A1(IJ,IJ,size) = A2(IJ,IJ,size+1);
%     [A1(IJ,IJ,size)] = preprocess(A1(IJ,IJ,size), 0);
    [A2(IJ,IJ,size), B2(1,IJ,size)] = OSR(A1(IJ,IJ,size));
    
    % также построить скорость нарастания энергии и приближение к исходной
    
%     A1(:,:,size) = A2(IJ,IJ,size+1);
%     A1(size+1:N1,:,size) = 0;
%     A1(:,size+1:N1,size) = 0;
%     [A2(:,:,size),B2(1,:,size)] = OSR(A1(:,:,size));
    
    a = A1(1,:,size+1);b = B2(1,:,size+1);
    for i = 2:N1
        a = [a A1(i,:,size+1)];
        b = [b B2(1,:,size+1)];
    end
    subplot(4,3,3*(N1-size)-2)
        plot(a),axis tight,grid on,set(gca,'xtick',[]),ylabel(num2str(size+1))
    subplot(4,3,3*(N1-size)-1)
        plot(B2(1,:,size+1)),axis tight,grid on,set(gca,'xtick',[])
    subplot(4,3,3*(N1-size)-0)
        plot(Eb2),axis tight,grid on,set(gca,'xtick',[])
end

c = zeros(1,N1);
for i = 1:N1
    c = c + B2(1,:,i);
end
figure,plot(1:N1,c,1:N1,B2(1,:,N1))

%% podporka
if length(IJ) == 2
    diag = [1/sqrt(2) 1/sqrt(2)];
    a = A1(IJ(1),IJ,i);
    proj(1) = a * diag';
    a = A1(IJ(2),IJ,i);
    proj(2) = a * diag';
    [~,ind] = min(proj);
    B(1,IJ(ind),i+2) = proj(ind);
    B(1,IJ,i+1) = proj - proj(ind);
    break
end

%% correlator
z = zeros(1,2*N1+1);
for k = 2:N1+1
    p = k-1;
    z(k) = 0;
    for j = 1:p
        z(k) = z(k) + S(j)*Et(j);
    end
end
for k = 2:N1-1
    z(N1+k) = 0;
    for j = k:N1
        z(N1+k) = z(N1+k) + S(j)*Et(j);
    end
end

figure
plot(z),axis tight,grid on

%% noise from preprocess
if nargin < 2, noise = 0; end  % by default or without 2nd arg
if noise == 0

else
    S = mean(a1);               % Среднее значение строк
    sh1 = 0.1*randn(N1);        % Шум

    e = S./sqrt(S*S');                          % Орт периодического сигнала
    
    sh2 = zeros(1,N1);
    for i = 1:N1
        sh2(i) = e*sh1(i,:)';                   % Проекция шума на орт
    end
    
    sh3 = zeros(1,N1);
    for i = 1:N1
        sh3(i,:) = e*sh2(i);                    % Коллинеарный шум
    end
    
    sh4 = zeros(1,N1);
    for i = 1:N1
        sh4(i,:) = sh1(i,:)-sh3(i,:);           % Ортогональный шум
    end
    
    a_2 = a1 + sh4;                             % Сигнал с ортогональным ему шумом
    a_3 = sh1;
    a_4 = sh3;
    a_5 = sh4;
    
end

%% was in preprocess
    for i = 1:N1
        a_ost1(i) = a1(i,:) * e';
        w1(i,:) = e * a_ost1(i);
        coll_ost1(i,:) = a1(i,:) - w1(i,:);
        perp_ost1(i,:) = w1(i,:) - S;
    end

    figure,plot(perp_ost','b'),hold,plot(coll_ost','.g'),plot(S,'r')

a_3 = coll_ost;
a_4 = perp_ost;
return

sum_ost_ort = zeros(1,N1);
vesa = zeros(1,N1);
for i = 1:N1
    sum_ost_ort(i) = mod_e + w(i);          % Сумма длин коллинеарных остатков и орта
    vesa(i) = sum_ost_ort(i) / mod_e;       % Отношение длин коллинеарных остатков и орта +1
end

norm_sum_ost_ort = zeros(N1);
norm_perp_ost = zeros(N1);
for i = 1:N1
    for j = 1:N1
        norm_sum_ost_ort(i,j) = ( S(j) + coll_ost(i,j) )/vesa(i);   % Нормировка коллинеарной составляющей
        norm_perp_ost(i,j) = perp_ost(i,j)/vesa(i);                 % Нормировка ортогональных остатков
    end
end

a_2 = zeros(N1);
for i = 1:N1
    a_2(i,:) = S + norm_perp_ost(i,:);
end


a_sum = zeros(N);
for i = 1:N
a_sum(i,:) = S + a_ort(i,:);
end

m_perp = mean(a_perp);

a_ort_dop = zeros(N1);
a_2 = zeros(N1);
for i = 1:N1
    a_ort_dop(i,:) = a_perp(i,:) - m_perp;
    a_2(i,:) = a_ort_dop(i,:);%c(i,:) + m_perp;
end

figure,plot(m_perp)
figure,plot(c')
figure,plot(a_ort_dop)

%% from orthogonalizations
figure(1)
subplot(311),plot(transform(A,'vector')),axis tight,grid on
subplot(312),plot(transform(b2,'vector_repeat')),axis tight,grid on
subplot(313),plot(transform(a2,'vector')),grid on,axis tight

figure(2)
subplot(511),plot([data(1,:) transform((axes(1,:)'*proj(1,:))','vector')]),grid on,axis tight
subplot(512),plot([data(2,:) transform((axes(2,:)'*proj(2,:))','vector')]),grid on,axis tight
subplot(513),plot([data(3,:) transform((axes(3,:)'*proj(3,:))','vector')]),grid on,axis tight
subplot(514),plot([data(4,:) transform((axes(4,:)'*proj(4,:))','vector')]),grid on,axis tight
subplot(515),plot([data(5,:) transform((axes(5,:)'*proj(5,:))','vector')]),grid on,axis tight

figure(3)
subplot(511),stem(other(1,:)),axis tight%,hold on, plot(sum(other,2),'.-k'),grid on
subplot(512),stem(other(2,:)),axis tight
subplot(513),stem(other(3,:)),axis tight
subplot(514),stem(other(4,:)),axis tight
subplot(515),stem(other(5,:)),axis tight

figure(4)
stem(len),grid on,axis tight

for i = 1:N
    for j = 1:N
        if i~=j
            SK(i,j) = data(i,:)*data(j,:)';
        end
    end
end
disp(max(max(abs(SK))))

for i = 1:N
    SK1(i) = b2*data(i,:)'; 
end
disp(max(abs(SK1)))


figure(1)
hold,subplot(311),plot(transform(A,'vector'),'r'),axis tight,grid on
hold,subplot(312),plot(transform(b2,'vector_repeat'),'r'),axis tight,grid on
hold,subplot(313),plot(transform(a2,'vector'),'r'),grid on,axis tight

figure(2)
hold on,subplot(511),plot([data(1,:) transform((axes(1,:)'*proj(1,:))','vector')],'r'),grid on,axis tight
hold on,subplot(512),plot([data(2,:) transform((axes(2,:)'*proj(2,:))','vector')],'r'),grid on,axis tight
hold on,subplot(513),plot([data(3,:) transform((axes(3,:)'*proj(3,:))','vector')],'r'),grid on,axis tight
hold on,subplot(514),plot([data(4,:) transform((axes(4,:)'*proj(4,:))','vector')],'r'),grid on,axis tight
hold on,subplot(515),plot([data(5,:) transform((axes(5,:)'*proj(5,:))','vector')],'r'),grid on,axis tight

figure(3)
hold on,subplot(511),stem(other(1,:),'r'),axis tight%,hold on, plot(sum(other,2),'.-k'),grid on
hold on,subplot(512),stem(other(2,:),'r'),axis tight
hold on,subplot(513),stem(other(3,:),'r'),axis tight
hold on,subplot(514),stem(other(4,:),'r'),axis tight
hold on,subplot(515),stem(other(5,:),'r'),axis tight

figure(4)
hold on,stem(len,'r'),grid on,axis tight

for i = 1:N
    for j = 1:N
        if i~=j
            SK(i,j) = data(i,:)*data(j,:)';
        end
    end
end
disp(max(max(abs(SK))))

%% surfaces from comparingABAxTri 13-11-14
figure(3)
subplot(131)
    surf(a.hist_a),axis tight
    title('a')
    xlabel('hist'),ylabel('noi')
subplot(132)
    surf(a.hist_ax),axis tight
    title('ax')
    xlabel('hist'),ylabel('noi')
subplot(133)
    surf(a.hist_tri),axis tight
    title('tri')
    xlabel('hist'),ylabel('noi')
figure(4)
subplot(131)
    surf(b.hist_a),axis tight
    title('a')
    xlabel('hist'),ylabel('noi')
subplot(132)
    surf(b.hist_ax),axis tight
    title('ax')
    xlabel('hist'),ylabel('noi')
subplot(133)
    surf(b.hist_tri),axis tight
    title('tri')
    xlabel('hist'),ylabel('noi')

figure
subplot(231)
    surf(P1_a),axis tight
    title('P1 a')
    xlabel('hist'),ylabel('comp')
subplot(232)
    surf(P1_ax),axis tight
    title('P1 ax')
    xlabel('hist'),ylabel('comp')
subplot(233)
    surf(P1_tri),axis tight
    title('P1 tri')
    xlabel('hist'),ylabel('comp')
subplot(234)
    surf(P2_a),axis tight
    title('P2 a')
    xlabel('hist'),ylabel('comp')
subplot(235)
    surf(P2_ax),axis tight
    title('P2 ax')
    xlabel('hist'),ylabel('comp')
subplot(236)
    surf(P2_tri),axis tight
    title('P2 tri')
    xlabel('hist'),ylabel('comp')

%%  from comparingABAxTri 13-11-18
% Histograms for different components
a.hist_b(inoi,:) = a.hist_b(inoi,:) + hist(a.st_cor_b(stat,:),ord);
a.hist_a(inoi,:) = a.hist_a(inoi,:) + hist(a.st_cor_a(stat,:),ord);
a.hist_ax(inoi,:) = a.hist_ax(inoi,:) + hist(a.st_cor_ax(stat,:),ord);
a.hist_tri(inoi,:) = a.hist_tri(inoi,:) + hist(a.st_cor_tri(stat,:),ord);

% Histograms
b.hist_b(inoi,:) = b.hist_b(inoi,:) + hist(b.st_cor_b(stat,:),ord);
b.hist_a(inoi,:) = b.hist_a(inoi,:) + hist(b.st_cor_a(stat,:),ord);
b.hist_ax(inoi,:) = b.hist_ax(inoi,:) + hist(b.st_cor_ax(stat,:),ord);
b.hist_tri(inoi,:) = b.hist_tri(inoi,:) + hist(b.st_cor_tri(stat,:),ord);

a.hist_b = a.hist_b/stn/N;
a.hist_a = a.hist_a/stn/N;
a.hist_ax = a.hist_ax/stn/N;
a.hist_tri = a.hist_tri/stn/N;

b.hist_b = b.hist_b/stn/N;
b.hist_a = b.hist_a/stn/N;
b.hist_ax = b.hist_ax/stn/N;
b.hist_tri = b.hist_tri/stn/N;

%% Old threshhold search
d = 1e-3;
for inoi = 1:noiN
   for i = 1:N
      
%       c(1,:) = abs( P1_a(inoi,i,:)-P2_a(inoi,i,:) );
%       [~,T_a(inoi,i)] = min( c(c>d) );
%       
%       c(1,:) = abs( P1_b(inoi,i,:)-P2_b(inoi,i,:) );
%       [~,T_b(inoi,i)] = min( c(c>d) );
      
      clear c
      c(1,:) = abs( P1_ax(inoi,i,:)-P2_ax(inoi,i,:) );
      [~,T_ax(inoi,i)] = min( c(c>d) );
      
%       if i==N,break;end
%       c(1,:) = abs( P1_tri(inoi,i,:)-P2_tri(inoi,i,:) );
%       [~,T_tri(inoi,i)] = min( c(c>d) );
      
   end
end

figure
plot(T_ax(:,1))

var1(1,:) = P1_ax(21,1,:);var2(1,:)=P2_ax(21,1,:);
figure,plot(var1,'x-b'),hold on,plot(var2,'+-g'),axis tight

%% between all in first window
b = length(Mov1{1}.period{1}.window{1}.ring);
c = length(Mov1{1}.period{1}.window{1}.ring{1}.per_comp(1,:));

img = zeros(b*c);
diag = img;
maxm = diag;
o = 0;
p = 0;
for j = 1:b
   for k = 1:c
      o = o+1;
      for m = 1:b
         for n = 1:c
            if p==b*c, p=0; end
            p = p+1;
            A = Mov1{1}.period{1}.window{i}.ring{j}.per_comp(k,:);
            B = Mov2{1}.period{1}.window{l}.ring{m}.per_comp(n,:);
            img(o,p) = (A/len(A)) * (B/len(B))';
            if o==p, diag(o,p) = img(o,p); end
            if 1-abs(img(o,p))<.1, maxm(o,p) = img(o,p); end
         end
      end
   end
end

% figure,surf(img)
% figure,surf(diag)
% figure,surf(maxm)
%% between continuous in half of windows
a = length(Mov1{1}.period{1}.window)/2;
b = length(Mov1{1}.period{1}.window{1}.ring);

o = 0;
p = 0;
q = 0;
img1 = zeros(a*b);
diag1 = zeros(1,a*b);
maxm1 = .5*ones(a*b);
for i = 1:a
   for j = 1:b
      o = o+1;
      for l = 1:a
         for m = 1:b
            if p==a*b, p=0; end
            p = p+1;
            A = Mov1{1}.period{1}.window{i}.ring{j}.per_comp(1,:);
            B = Mov2{1}.period{1}.window{l}.ring{m}.per_comp(1,:);
            img1(o,p) = (A/len(A)) * (B/len(B))';
            if o==p, q=q+1; diag1(q) = img1(o,p); end
            if 1-abs(img1(o,p))<.1, maxm1(o,p) = img1(o,p); end
         end
      end
   end
end

img1 = img1/2 + .5;
figure,imshow(img1)
figure,stem(diag1),axis tight
maxm1 = maxm1/2 + .5;
figure,imshow(maxm1)
%% between same continuous in all windows
a = length(Mov1{1}.period{1}.window);
b = length(Mov1{1}.period{1}.window{1}.ring);

o = 0;
cor = zeros(1,a*b);
rinI = cor;
for i = 1:a
   for j = 1:b
      o = o+1;
      A = Mov1{1}.period{1}.window{i}.ring{j}.per_comp(1,:);
      B = Mov2{1}.period{1}.window{i}.ring{j}.per_comp(1,:);
      cor(o) = (A/len(A)) * (B/len(B))';
   end
   rinI(o) = 1;
end

figure,stem(rinI,'.:r')
hold on,stem(abs(cor),'.k'),axis tight
