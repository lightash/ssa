% close all
figure
J = 36;
winL2 = winL/2;
Blen = zeros(J,winL2);
A1 = cell(J,winL2);
A2 = cell(J,winL2);
B2 = cell(J,winL2);
for j = 1:J
   for i = 1:winL2
      diap = winL2*(j-1)+1:i+(j-1)*winL2;
      A1{j,i} = nrm(f(Bnum{1}(diap),winL2+1:winL2+i));
      [B2{j,i},A2{j,i}] = OSR(A1{j,i});
%       Bdim(i) = log10(max(abs(B2{j,i})));
      Blen(j,i) = len(B2{j,i});
   %    figure(1),plot(B2{i},'.-'),hold on
   %    figure(2),plot(A1{i}','.-'),hold on
   end
   plot(Blen(j,:),'.-'),hold on
end
%%
plot(mean(Blen,1),'.-r'),hold on
%%
plot(1./(1:winL2),'k')