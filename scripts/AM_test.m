clc;
close all;
clear all;

load('d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\S001\R03\S001R03')

T = floor(sqrt(size(mov1{1},2)));
F3 = 32;

signal = transform(mov1{F3}(1,1:T^2),'matrix');

[e,c,a2] = AM(signal);

new = signal - ( e' * (c-mean(c)) )';

figure
plot(transform(signal,'vector'),'g'),hold on;
plot(transform( (e'*c)' ,'vector'),'r'),hold on;
plot(transform(e*mean(c),'vector_repeat'),'b'),hold on;
plot(transform(new,'vector'),'k'),hold on;
plot(transform(a2,'vector'),'y'),hold on;
legend('signal','(e''*c)''','e*mean(c)','new','a2')

[e1,c1,a21] = AM(new);

figure
for i=1:T
   plot([1:T]+T*(i-1),new(i,:),'g'),hold on;
	plot([1:T]+T*(i-1),e1*c1(i),'.r'),hold on;
	plot([1:T]+T*(i-1),e1*mean(c1),'b'),hold on;
% 	plot([1:T]+T*(i-1),a2(i,:),'g'),hold on;
end

figure,plot(a21'-a2')

% for i = 1:T
%    orth(i) = log10(abs( (e*c(i)) * a2(i,:)' ));
% end
% disp(['Orthogonality is more than 10^' num2str(max(orth))])
% plot(sum(a2,1))