clc;
close all;
clear all;

% ��������
img = imread('Images\img2.png');
grayimg = rgb2gray(img);
N = 256;
n = 1;
a = double(grayimg( 1: n: n*N, 1: n: n*N ));


[sr_mat, a_coll, ~, l_coll] = preprocess(a);


figure('Color','w')
subplot(121),imshow(uint8(255*pic_norm(a,'to_one')))
    title({'��������';'�����������'},'FontName','Times New Roman','FontSize',14)

figure('Color','w')
subplot(121),imshow(uint8(255*pic_norm(sr_mat,'to_one')))
    title('|| �_�_� ||','FontName','Times New Roman','FontSize',14)
subplot(122),plot(255*pic_norm(l_coll,'to_one'),'k'),axis tight,grid on
    title('�_�_�_�_�','FontName','Times New Roman','FontSize',14)

figure('Color','w')
subplot(121),imshow(uint8(255*pic_norm(sr_mat+a_coll,'to_one')))
    title('|| �_�_� || + || �_�_�_�_� ||','FontName','Times New Roman','FontSize',14)
