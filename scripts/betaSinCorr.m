clc;
close all;
clear all;

pat = 1;
load(['d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\S00' num2str(pat) '\R03\S00' num2str(pat) 'R03']) % R03|07|11
F = 1:.1:80;%[14:18 (10:14)*4/3 (19:26)*5/7];
Ts = 1./F;
meancoef = 4;
meannum = 16;
winS = 12;
Nphas = 360/15;
WIN = 40;
%       inRaw = raw{32};
for freq = 1:length(F)
% figure
   
   t = linspace(0,winS/Fd,winS); % 0 : 1/Fd : Ts(freq);
   
   for movm = 1
      disp([num2str(F(freq)) '   ' num2str(movm-1)])

      eval(['inRaw = mov' num2str(movm-1) '{32};'])
%       inRaw = raw{32};
      
      Cor = zeros(size(inRaw));
      defi = zeros(size(inRaw));
      M = zeros(size(inRaw));
      winI = WIN-winS/2:WIN+winS/2-1;
      for win = winI %1:size(inRaw,2) - length(t)+1

         Co = zeros(1,size(inRaw,1));
         dfi = zeros(1,size(inRaw,1));
         for rept = 1:size(inRaw,1)

            inwin = inRaw(rept, win:win + length(t)-1 );

            M(rept,win ) = mean(inwin);%+ winS/2

            S = zeros(Nphas,length(t));
            C = zeros(1,Nphas);
            for phas = 1:Nphas
               fi(phas) = 2*pi*(phas-1)/Nphas;
               S(phas,:) = (sin( 2*pi*F(freq)*t + fi(phas) ));
               C(phas) = (inwin - M(win )) * S(phas,:)';%+ winS/2
            end
            
            [Co(rept), phase_index] = max(C);
            dfi(rept) = fi(phase_index);
         end

         Cor(:,win ) = Co';%+ winS/2
         defi(:,win ) = dfi';%+ winS/2
      end
      
      Corr{freq,movm}(:,WIN) = mean(Cor(:,winI),2);
      Deltafi{freq,movm} = defi;
      Meanwin{freq,movm} = M;
      
      % Recovery
%       for rept = 1:size(Corr{freq,movm},1)
%          Rec{freq,movm}(rept,:) = Corr{freq,movm}(rept,:) .* sin(2*pi*F(freq)*t(winS/2+1) + Deltafi{freq,movm}(rept,:))...
%             + Meanwin{freq,movm}(rept,:);
%       end
%       plot(linspace(0,length(Rec{freq,movm})/Fd,length(Rec{freq,movm})),Rec{freq,movm}(rept,:),'r'),axis tight,ylim([-170 170])

   end
   
%    inRaw = inRaw - Rec{freq,movm};
end

figure
for freq = 1:length(F)
   plot(F(freq),Corr{freq,movm}(:,WIN )-300*(0:14)','k'),hold on%+ winS/2
%    subplot(length(F),1,freq),plot(inRaw,'g')
%    hold on,plot(Corr{freq,1}+Meanwin{freq,1}),ylim([-200 400]),ylabel(F(freq))
end

% % Vertical lines between moves
% mark = [0 0 -150 150];
% for ann = 1:length(annot)
%    if annot(ann) > 0
%       mark = mark + [length(mov0{1})/Fd length(mov0{1})/Fd 0 0];
%    else
%       mark = mark + [length(mov1{1})/Fd length(mov1{1})/Fd 0 0];
%    end
%    plot(mark(1:2),mark(3:4),'.-k')
% end









