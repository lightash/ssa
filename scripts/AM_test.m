clc;
close all;
clear all;

load('d:\Dropbox\Signals\Msasr.mat')

M = size(Msasr,1);
N = size(Msasr,2);

[e,c,a2] = AM(Msasr);

figure,for i=1:26,plot([1:1698]+1698*(i-1),Msasr(1,:),'b'),hold on;end
for i=1:26,plot([1:1698]+1698*(i-1),e*c(i),'r'),hold on;end
for i=1:26,plot([1:1698]+1698*(i-1),a2(i,:),'g'),hold on;end

figure,plot(a2')

for i = 1:M
   orth(i) = log10(abs( (e*c(i)) * a2(i,:)' ));
end
disp(['Orthogonality is more than 10^' num2str(max(orth))])
plot(sum(a2,1))