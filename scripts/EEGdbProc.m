clc
close all
clear all

load('d:\Dropbox\MATLAB\������\Signals\EEG Motor Movement-Imagery Dataset\Processed\S001R03')

T = 2:floor(sqrt(min(diff(mar))));
F3 = 32;

% Movement 1 portraits
N = [25 24 23];
for sig = 1:size(mov1{F3},1)
   signal = sil{F3}(sig,:);
   SilPort{sig} = perPort(signal,N);
   
   signal = mov1{F3}(sig,:);
   Mov1Port{sig} = perPort(signal,N);
end
disp('Done')

%% informativity in one patient between same continuous for all movements in all windows
a = length(Mov1Port);
b = length(Mov1Port{1}.period{1}.window);
c = length(Mov1Port{1}.period{1}.window{1}.svproj(1,:));

st1 = zeros(c,a*b*c);
st2 = st1;
for ci = 1:c
   i = 1;
   for ai = 1:a
      for bi = 1:b
         st1(ci,i:i-1+c) = SilPort{ai}.period{1}.window{bi}.svproj(:,ci)';
         st2(ci,i:i-1+c) = Mov1Port{ai}.period{1}.window{bi}.svproj(:,ci)';
%          if sum(isnan(st1(ci,i:i-1+c)))>0, disp([0 ai bi ci]); end
         i = i+c;
      end
   end
   disp([ci c])
end

close all
for i=1:N
   figure(i)
   subplot(311),stem(st1(i,:),'.')
   subplot(312),stem(st2(i,:),'.')
end


ord = -50:1:50;
for i=1:N
   figure(i)
   subplot(313),plot(ord,hist(st1(i,:),ord),ord,hist(st2(i,:),ord),'r')
end








