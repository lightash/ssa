clc;
close all;
clear all;

load('d:\Dropbox\MATLAB\Био. комп. сис\signal\Ms.mat')
% load('D:\Dropbox\MATLAB\Диплом\Signals\Raw\EEG_from_Popov\healthy_background\DX_EEG_BEAfon_AF_free.mat')
%     for z=1:31,noi = 0:.1:3;disp(z)
N = 64;%32;
n = 16;%32;
A = Ms(1:N,1:n:n*N);%+noi(z)*randn(64);
% A=transform(h.data(1,1:N*N),'matrix');
% [Bimp,lowTr] = imp_OSR(A,'from_end');
% Bimp = Bimp';

[a2,b2] = OSR(A);

data = a2;%Bimp;

for axi = 1:N
%     disp(axi)
    
    % Normalization
    axes(axi,:) = data(axi,:);
    len = sqrt( axes(axi,:)*axes(axi,:)' );
    axes(axi,:) = axes(axi,:) / len;

    for i = axi+1:N
        proj(axi,i) = data(i,:) * axes(axi,:)';  % Lenghts of projections on chosen axes
        vproj(i,:) = proj(axi,i) * axes(axi,:);  % Vectors of projections
        data(i,:) = data(i,:) - vproj(i,:);
        data(axi,:) = data(axi,:) + vproj(i,:);
    end
    
end

% figure('Color','w')
% subplot(311),plot(A'),axis tight,grid on
% subplot(312),plot(a2'),grid on,axis tight
% subplot(313),plot(data'),grid on,axis tight

for i = 1:N
    for j = 1:N
        if i~=j
            SK(i,j) = data(i,:)*data(j,:)';
        end
    end
end
disp(max(max(abs(SK))))

% for i=1:64,pow2(z,i) = data(i,:)*data(i,:)';pow1(z,i) = Bimp(i,:)*Bimp(i,:)';end

% figure,plot(pow1(z,:)),hold on,plot(pow2(z,:),'g'),axis tight
% legend(num2str(sum(pow1(z,:))), num2str(sum(pow2(z,:))))

%     end
    
% for i=1:z
% p1(i)=sqrt(sum( (pow1(1,:)-pow1(i,:)).^2 ))/sum(pow1(1,:));
% 
% p2(i)=sqrt(sum((pow2(1,:)-pow2(i,:)).^2))/sum(pow2(1,:));
% end
% 
% figure,plot(p1,'.-b')
% hold on,plot(p2,'.-g'),axis tight
% legend('Bimp','data')

for i=1:N, lo(i)=sqrt( data(i,:)*data(i,:)' ); end
for i=1:N, loun(i)=sqrt( a2(i,:)*a2(i,:)' ); end
figure(3),stem(lo,'b')
a1 = zeros(N);
for i=1:N
    ao(i,:) = a2(i,:)/sqrt( a2(i,:)*a2(i,:)' );
    for j=1:N
        a1(i,:)=a1(i,:)+( data(j,:) * ao(i,:)' ) * ao(i,:);
    end
end
% figure,plot(( a1/( max(max(a1))/max(max(a2)) ) )'),axis tight
% figure,plot(a2'),axis tight
    for z=1:16
A = Ms(1:N,1:n:n*N)+.5*randn(64);
% A=transform(h.data(1,1:N*N),'matrix')+.5*randn(64);

[a2,b2] = OSR(A);
data = a2;
for i=1:N, lun(z,i)=sqrt( a2(i,:)*a2(i,:)' ); end
% figure('Color','w')
% subplot(311),plot(A'),axis tight,grid on
% subplot(312),plot(a2'),grid on,axis tight

for axi = 1:N
    axes(axi,:) = data(axi,:);
    len = sqrt( axes(axi,:)*axes(axi,:)' );
    axes(axi,:) = axes(axi,:) / len;

    for i = axi+1:N
        proj(axi,i) = data(i,:) * axes(axi,:)';  % Lenghts of projections on chosen axes
        vproj(i,:) = proj(axi,i) * axes(axi,:);  % Vectors of projections
        data(i,:) = data(i,:) - vproj(i,:);
        data(axi,:) = data(axi,:) + vproj(i,:);
    end
end
% subplot(313),plot(data'),grid on,axis tight

for i=1:N, l(z,i)=sqrt( data(i,:)*data(i,:)' ); end
    end
figure(3),hold on,stem(mean(l,1),'r')
hold on,stem((mean(l,1)-lo),'g'),axis tight
sum(abs((mean(l,1)-lo)))
sum(abs((mean(lun,1)-loun)))