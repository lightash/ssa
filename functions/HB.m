function [I] = HB(p, q)
% Histograms border informativity criterium.
%    p and q are (num,len) matrices of data;
%    I is 1 by num_hists vector of informativity numbers.

num_hists = size(p,2);

Nbars = 16;
for i = 1:num_hists
   [~,scale(i,:)] = hist( [p(:,i);q(:,i)] ,Nbars);
   
   Mp(i) = mean(p(:,i));
   Hp(:,i) = hist( p(:,i) ,scale(i,:))/size(p,1)';
   % Распределить наименьший столбец по нулям
   [minh,indh] = min(Hp( Hp(:,i)>0, i));
   mzh = minh/( length( Hp( Hp(:,i)==0, i)) +1);
   Hp( Hp(:,i)==0 | Hp(:,i)==minh, i) = mzh;
   Hp(:,i) = Hp(:,i)/sum(Hp(:,i));
   
   Mq(i) = mean(q(:,i));
   Hq(:,i) = hist( q(:,i) ,scale(i,:))/size(q,1)';
   % Распределить наименьший столбец по нулям
   [minh,indh] = min(Hq( Hq(:,i)>0, i));
   mzh = minh/( length( Hq( Hq(:,i)==0, i)) +1);
   Hq( Hq(:,i)==0 | Hq(:,i)==minh, i) = mzh;
   Hq(:,i) = Hq(:,i)/sum(Hq(:,i));
end

Nbars = 256;
border = zeros(1,Nbars);
for i = 1:num_hists
   leftN(i) = Mp(i) < Mp(i);
   
   scali(i,:) = linspace(scale(i,1),scale(i,end),Nbars);
   
   Hip(:,i) = interp1(scale(i,:),Hp(:,i),scali(i,:));
   Hip(:,i) = Hip(:,i)/sum(Hip(:,i));
   Hiq(:,i) = interp1(scale(i,:),Hq(:,i),scali(i,:));
   Hiq(:,i) = Hiq(:,i)/sum(Hiq(:,i));
   
   for bari = 1:Nbars
      if leftN(i)
         brd(bari) = sum(Hip( 1:bari ,i)) - sum(Hiq( bari+1:end ,i));
      else
         brd(bari) = sum(Hiq( 1:bari ,i)) - sum(Hip( bari+1:end ,i));
      end
   end
   [~,border(i)] = min(abs(brd));
   bord_pq(i) = scali(i,border(i));
   
   if leftN(i)
      p00(i) = sum(Hip(1:border(i),i));
      p11(i) = sum(Hiq(border(i)+1:end,i));
   else
      p00(i) = sum(Hip(border(i)+1:end,i));
      p11(i) = sum(Hiq(1:border(i),i));
   end
end
I = mean([p00;p11],1);