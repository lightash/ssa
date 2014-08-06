% function [frot] = ax_rot(f,sk,ss)
clc;
close all;
clear all;

load('Signals\Msasr')
f = Msasr(1,:);

% f = exp(0:-1e-1/4:-10/4)-exp(0:-1e-1/2:-10/2);

sk = .05;       % Шаг квантования
ss = 1.0;       % Шаг сетки (коэф. сжатия)

x = 1:length(f);
xi = 1:ss:length(f);
fi = interp1(x,f,xi);

xk = (round(min(fi)/sk)-1)*sk : sk : (round(max(fi)/sk)+1)*sk;

l = 0;
for i = 1:length(fi)
    k = abs( fi(i)-xk );
    if min(k) < ss/2
        l = l+1;
        fk(l) = sk*round(fi(i)/sk);
        xn(l) = i;
    end
end

fk2(1) = fk(1);
j=1;
for i = 2:length(fk)
    if fk(i)-fk(i-1) ~= 0
        j=j+1;
        fk2(j) = fk(i);
        xk2(j) = i;
    end
end
% figure,plot(xk2,'.')
% figure,plot(1:20,xk2(1:20),'.')

frot(1) = xk2(1);
loc_ex = 0;
for i = 2:length(fk2)
    if fk2(i) > fk2(i-1)
        
        frot(i) = loc_ex+xk2(i);
        loc_ex = frot(i);
        
    elseif fk2(i) < fk2(i-1)
        
        frot(i) = loc_ex-xk2(i);
        loc_ex = frot(i);
        
    else
        frot(i) = frot(i-1);
    end
end

figure
subplot(211),plot(xi,fk,'.',x,f,'r'),grid,axis tight
subplot(212),plot(xk2,frot,'.-'),grid,axis tight
