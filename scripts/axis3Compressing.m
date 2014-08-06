clc;
close all;
clear all;

% ��������
img = imread('Images\img0.png');
% grayimg = img(:,:,1);
grayimg = rgb2gray(img);
N1 = 128;%272;%350;
n = 1;
a = double(grayimg(1:n:n*N1,1:n:n*N1));
% a = pic_norm(a);
% a = preprocess(a,0);

% % % ������������ �������
% % N1 = 4;
% % a = [1 6 2 4 ; 3 1 5 1 ; 6 2 2 3 ; 1 3 2 4];

figure(1)%('Color','w')
subplot(241),imshow(uint8(pic_norm(a))),title('a')

% % pic = a(1,:);
% % for i = 2:N1
% %     pic = [pic a(i,:)];
% % end
% % subplot(231),plot(pic),grid on,axis tight,title('a')

[~,b] = OSR(a,1);
bo = b/sqrt(b*b');

a_cp = mean(a);

c = a_cp - b;
co = c/sqrt(c*c');

di = zeros(N1);
for i = 1:N1
    di(i,:) = a(i,:) - b - c;
end

di_c = ( di(1,:) * co' ) * co;
di_d = di(1,:) - di_c;
d = di_d / sqrt(di_d*di_d');
do = d/sqrt(d*d');

figure(2)
subplot(241),imshow(uint8(pic_norm(di))),title('di')

cm = transform(c,'matrix_repeat');
dm = transform(d,'matrix_repeat');
bm = transform(b,'matrix_repeat');
a_cpm = transform(a_cp,'matrix_repeat');

figure(1)
subplot(242),imshow(uint8(pic_norm(cm))),title('c')
subplot(243),imshow(uint8(pic_norm(dm))),title('d')
subplot(244),imshow(uint8(pic_norm(bm))),title('b')


ac = zeros(1,N1);ad = ac;ab = ad;
for i = 1:N1
    ac(i) = a(i,:)*co';
    ad(i) = a(i,:)*do';
    ab(i) = a(i,:)*bo';
end

acm = transform(ac,'matrix_repeat');
adm = transform(ad,'matrix_repeat');
abm = transform(ab,'matrix_repeat');

figure(2)
subplot(242),imshow(uint8(pic_norm(acm))),title('ac')
subplot(243),imshow(uint8(pic_norm(adm))),title('ad')
subplot(244),imshow(uint8(pic_norm(abm))),title('ab')

cc = c*co';
cd = c*do';
cb = c*bo';
dc = d*co';
dd = d*do';
db = d*bo';
bc = b*co';
bd = b*do';
bb = b*bo';
a_cpc = a_cp*co';
a_cpd = a_cp*do';
a_cpb = a_cp*bo';
dic = zeros(1,N1);did = dic;dib = did;
for i = 1:N1
    dic(i) = di(i,:)*co';%ac(i) - bc - cc;
    did(i) = di(i,:)*do';%ad(i) - bd - cd;
    dib(i) = di(i,:)*bo';%ab(i) - bb - cb;
end

% Return to N1-dim

acn = zeros(N1);adn = acn;abn = acn;
for i = 1:N1
    acn(i,:) = ac(i)*co;
    adn(i,:) = ad(i)*do;
    abn(i,:) = ab(i)*bo;
end

figure(2)
subplot(246),imshow(uint8(pic_norm(acn))),title('acn')
subplot(247),imshow(uint8(pic_norm(adn))),title('adn')
subplot(248),imshow(uint8(pic_norm(abn))),title('abn')

an = zeros(N1);
for i = 1:N1
    an(i,:) = acn(i,:) + adn(i,:) + abn(i,:);
end


a3c = dic + sqrt(c*c');
a3d = did + sqrt(d*d');
a3b = dib + sqrt(b*b');

a3n = zeros(N1);
for i = 1:N1
    a3n(i,:) = a3c(i)*co' + a3d(i)*do' + a3b(i)*bo';
end

figure(3),imshow(uint8(pic_norm(a3n)))


figure(1)
subplot(245),imshow(uint8(pic_norm(an))),title('an')

cn = zeros(N1);bn = cn;dn = cn;din = cn;
for i = 1:N1
    cn(i,:) = cc*co + cd*do + cb*bo;
    dn(i,:) = dc*co + dd*do + db*bo;
    bn(i,:) = bc*co + bd*do + bb*bo;
    din(i,:) = dic(i)*co + did(i)*do + dib(i)*bo;
end

figure(1)
subplot(246),imshow(uint8(pic_norm(cn))),title('cn')
subplot(247),imshow(uint8(pic_norm(dn))),title('dn')
subplot(248),imshow(uint8(pic_norm(bn))),title('bn')

figure(2)
subplot(245),imshow(uint8(pic_norm(din))),title('din')

% figure
% subplot(121),imshow(uint8(pic_norm(a)*255)),title('a')
% subplot(122),imshow(uint8(pic_norm(an)*255)),title('an')








