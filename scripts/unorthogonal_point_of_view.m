clc;

btypeN = 2;
perN = 2*length(ann{2});
fg = [f(anniNAV{1}(1:perN/2),:);f(anniNAV{2}(1:perN/2),:)];

N = port{1};
A = port{2};

eN = nrm(N - A);
eA = nrm(A - N);
nor = nrm(A - (A - N)/2);

% unor{1} = eN * (eN * N')';
% unor{2} = eA * (eA * A')';


kN = 1+120;
desN = zeros(kN,btypeN);
desA = zeros(kN,btypeN);
for k = 1:kN
   phi(k) = (k-1)*2*pi/(kN-1);
   
   unor{2} = eA*cos(phi(k)) + nor*sin(phi(k));
   unor{1} = -unor{2};
   
   for per = 1:perN
      sig = fg(per,:);
      sig = nrm(sig - mean(sig));

%       fgN(per,:) = unor{1} * (sig * unor{1}')';
%       fgA(per,:) = unor{2} * (sig * unor{2}')';
      for btype = 1:btypeN
         cor(btype,per) = sig * unor{btype}';
   %       cor(btype,per) = (cor(btype,per) +1)/2;
      end

      [~,ind] = max(cor(:,per));
      if per <= length(ann{2})
         desN(k,ind) = desN(k,ind) + 2/perN;
      else
         desA(k,ind) = desA(k,ind) + 2/perN;
      end
   end
end
% hold on,plot(unor{1})
% hold on,plot(unor{2},'g')
% hold on,plot(sig,'r')
%%
% figure
% plot(subspN,'-'),hold on
% plot(subspA,'-g'),hold on
% plot(unorN,':'),hold on
% plot(unorA,':g'),hold on

figure
subplot(btypeN,btypeN,1),plot(phi/pi,desN(:,1),'.-'),axis tight,ylim([0 1])
subplot(btypeN,btypeN,2),plot(phi/pi,desN(:,2),'.-'),axis tight,ylim([0 1])
subplot(btypeN,btypeN,3),plot(phi/pi,desA(:,1),'.-'),axis tight,ylim([0 1])
subplot(btypeN,btypeN,4),plot(phi/pi,desA(:,2),'.-'),axis tight,ylim([0 1])


















