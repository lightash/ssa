clc;
close all;
clear all;
% tic
load('Signals\Msasr')
n = 50;
input = [transform(Msasr(:,1:n:n*26),'vector') Msasr(1,1:n:n*16)];

N = floor(sqrt(length(input)));  % Chosen period %26
N1 = N^2;
Shifts = length(input) - N1;  % Number of shifts %16
Num_iter = 4;  % Number of iterations

A3 = zeros(N,N,Shifts+1,Num_iter);
Init = zeros(Num_iter,length(input));
windows = zeros(Shifts+1,N,Num_iter);
win_en = zeros(Shifts+1,Num_iter);
Rep = zeros(Num_iter,N1);

for per = N:-1:19
    disp(['period = ' num2str(per)])
    
    data = input;
    for iter = 1:Num_iter  % Just repeated iterations
        disp(['    iteration = ' num2str(iter)])

        Init(iter,:) = data;  % Initial data for iteration

        [windows(:,:,iter)] = win_OSR(data,per);

        Rep(iter,:) = transform(mean(windows(:,:,iter),1),'vector_repeat');  % Obtained repeated components
        dd = [Rep(iter,:) Rep( iter,1:length(data)-length( Rep(iter,:) ) )];
        data = data - dd;
    end

end

figure
plot(sum(Rep,1),'r'),hold on,plot(Init(1,:),'b'),axis tight,grid on
% toc

