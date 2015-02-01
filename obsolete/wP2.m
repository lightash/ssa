clc;
% close all;
clear all;

load('Signals\Msasr')
N = 26;
s = 1;%503;
n = 50;%1;
a1 = Msasr(1:N,s:n:s+n*N-1);
N1 = N*N;

ss = transform(a1,'vector');

z = zeros(N+1,2*N-1,11,20);
z_imp = zeros(N+1,2*( 2*N-1 )-1,N,11,20);
for iter = 1:20
    disp(['noise = ' num2str(iter)])
for k = 1:11
    disp(['    noise_ampl = ' num2str(.2*(k-1))])
    shum = .2*(k-1)*randn(1,N1);
    S = ss + shum;
    
    [B, A3] = imp_OSR(a1);
    
%     for imp = 1:N
%         for i = 1:N
%             Et_imp(imp,:) = xcorr(sum(B(1,:,1:imp),3),ss( (i-1)*N+1 : i*N ));
%             ENEt_imp = Et_imp(imp,:)*Et_imp(imp,:)';
%             Et_imp(imp,:) = Et_imp(imp,:)/sqrt(ENEt_imp);
%             
%             sig(i,:) = xcorr(S( (i-1)*N+1 : i*N ),sum(B(1,:,1:imp),3));
%             ENsig = sig(i,:)*sig(i,:)';
%             sig(i,:) = sig(i,:)/sqrt(ENsig);
%             
%             z_imp(i,:,imp,k) = xcorr(Et_imp(imp,:),sig(i,:));
%             z_imp(N+1,:,imp,k,iter) = z_imp(N+1,:,imp,k,iter) + z_imp(i,:,imp,k);
%         end
% 
%         fmax_imp(k,imp,iter) = z_imp(N+1,26,imp,k,iter);
%         smax_imp(k,imp,iter) = max([z_imp(N+1,1:25,imp,k,iter) z_imp(N+1,27:51,imp,k,iter)]);
%     end

    Et = sum(B(1,:,:),3);
    ENEt = Et*Et';
    Et = Et/sqrt(ENEt);

    for i = 1:N
        sig(i,:) = S( (i-1)*N+1 : i*N );
        ENsig = sig(i,:)*sig(i,:)';
        sig(i,:) = sig(i,:)/sqrt(ENsig);
        z(i,:,k) = xcorr(Et,sig(i,:));
        z(N+1,:,k,iter) = z(N+1,:,k,iter) + z(i,:,k);
    end

    fmax(k,iter) = z(N+1,26,k,iter);
    smax(k,iter) = max([z(N+1,1:25,k,iter) z(N+1,27:51,k,iter)]);
end
end
figure
subplot(411),plot(S),axis tight,grid on
subplot(412),plot(transform(Et,'vector_repeat')),axis tight,grid on
% subplot(413),plot(z(1:N,:)'),axis tight,grid on
subplot(413),plot(mean(fmax,2),'b'),hold on,plot(smax,'r'),axis tight,grid on
for j = 1:size(z,3), zz(:,j) = mean(z(N+1,:,j,iter),4); end
subplot(414),plot(zz),axis tight,grid on

% figure
% for i = 1:6
%     for j = 1:size(z_imp,4), zz(:,j) = mean(z_imp(N+1,:,i,j,iter),5); end
%     subplot(6,2,2*i-1),plot(zz),axis tight,grid on
%     subplot(6,2,2*i),plot(mean(fmax_imp(:,i,:),3),'b'),hold on,plot(mean(smax_imp(:,i,:),3),'r'),axis tight,grid on
% end
% figure
% for i = 1:6
%     for j = 1:size(z_imp,4), zz(:,j) = mean(z_imp(N+1,:,i+20,j,iter),5); end
%     subplot(6,2,2*i-1),plot(zz),axis tight,grid on
%     subplot(6,2,2*i),plot(mean(fmax_imp(:,i+20,:),3),'b'),hold on,plot(mean(smax_imp(:,i+20,:),3),'r'),axis tight,grid on
% end

