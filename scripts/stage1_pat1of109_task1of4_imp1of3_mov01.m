clc;
% close all;
clear all;

load('d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\S001\R03\S001R03')

T = floor(sqrt(min(size(mov1{1},2))));
F3 = 32;

signal = mov1{F3}(1,1:T^2);

n = 400;
remPort = cell(1,n);
a1 = signal;
e = zeros(1,n);
for i = 1:n
   disp(i/n)
   remPort{i} = perPort(a1,T,'AM');
   for j = 1:T
      or = remPort{i}{1}{1}.svproj / len(remPort{i}{1}{1}.svproj);
      b2(j,:) = or * (or * remPort{i}{1}{1}.rem(j,:)');
      e(i) = e(i) + len(b2(j,:))/T;
   end
   a1 = a1 - transform(b2,'vector');
end
%%
figure,plot(e(1:n)),axis tight
% figure
% for i = 1:n
%    subplot(n,1,i),plot(remPort{i}{1}{1}.svproj,'.-'),axis tight
% end


