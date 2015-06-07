clc;
close all;
clear all;

load('d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\S001\R03\S001R03')

T = floor(sqrt(min(size(mov1{1},2))));
F3 = 32;
Nmov0 = size(mov0{F3},1);
Nmov1 = size(mov1{F3},1);
Nmov2 = size(mov2{F3},1);
for mov = 1:Nmov0+Nmov1+Nmov2
   if mov <= Nmov0
      signal = mov0{F3}(mov,1:T^2);
   elseif mov <= Nmov0 + Nmov1
      signal = mov1{F3}(mov-Nmov0,1:T^2);
   elseif mov <= Nmov0 + Nmov1 + Nmov2
      signal = mov2{F3}(mov-Nmov0-Nmov1,1:T^2);
   end
   sigmat = transform(signal,'matrix');

   a1 = sigmat;
%    n = 1;
%    e = zeros(1,n);
%    for iter = 1:n
%       disp([num2str(mov) ' ' num2str(iter/n,2)])
      [E,C,A3{mov}] = impAM(a1);
      or(mov,:) = sum( GSOrth(E) ,1);
%       e(iter) = len(or(mov,:));
%       a1 = A3;
%       for imp = 1:T
%          a1 = a1 + ( C(imp,:) )' * ( E(imp,:)/len(E(imp,:)) );
%       end
%    end
end

Lmov0 = length(mov0{F3}(1,:));
Lmov1 = length(mov1{F3}(1,:));
signal = [];sigor = [];siga3 = [];
i0=0;i1=0;i2=0;
for ann = 1:length(annot)
   if annot(ann) == 0
      i0 = i0 + 1;
      signal = [signal mov0{F3}(i0,:)];
      sigor = [sigor transform(or(i0,:),'vector_repeat') zeros(1,Lmov0-T^2)];
      siga3 = [siga3 transform(A3{i0},'vector') zeros(1,Lmov0-T^2)];
   elseif annot(ann) == 1
      i1 = i1 + 1;
      signal = [signal mov1{F3}(i1,:)];
      sigor = [sigor transform(or(Nmov0+i1,:),'vector_repeat') zeros(1,Lmov1-T^2)];
      siga3 = [siga3 transform(A3{i0},'vector') zeros(1,Lmov1-T^2)];
   elseif annot(ann) == 2
      i2 = i2 + 1;
      signal = [signal mov2{F3}(i2,:)];
      sigor = [sigor transform(or(Nmov0+Nmov1+i2,:),'vector_repeat') zeros(1,Lmov1-T^2)];
      siga3 = [siga3 transform(A3{i0},'vector') zeros(1,Lmov1-T^2)];
   end
end

figure
plot(signal),hold on
plot(sigor,'g'),hold on
plot(siga3,'r'),hold on
plot(sigor+siga3,'k'),hold on
axis tight


















