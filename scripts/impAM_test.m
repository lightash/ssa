clc;
close all;
clear all;

load('d:\Dropbox\Signals\Msasr.mat')

M = size(Msasr,1);
N = M;%size(Msasr,2);

n = 30;
in = Msasr(:,1:n:N*n);
[E,C,A3] = impAM(in);

figure
for i = 1:M
   plot([1:N]+N*(i-1),in(1,:),'b'),hold on;
   for j = 1:N
      plot([1:N]+N*(i-1),E(i,:)*C(i,j),'r'),hold on;
%    plot([1:1698]+1698*(i-1),a2(i,:),'g'),hold on;
   end
end

% figure,plot(A3')