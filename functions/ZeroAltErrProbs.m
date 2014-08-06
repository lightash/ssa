function [a, b, c] = ZeroAltErrProbs(a1, b1, noi, ord)
%ABAxTriDecomp   Calculating probabilities of
% zero and alternative hypothesis and
% errors of the first and second kind
%   a1 - structure with first object variables
%   b1 - structure with second object variables
%   noi - noise vector
%   ord - positions of histogram bars

a = a1;
b = b1;

N = size(a.a,1);
noiN = length(noi);
ordN = length(ord);

% 1) Нормировка площади плотности вероятности
Norm1 = zeros(noiN,N,8);
for inoi = 1:noiN
   for i = 1:N
      
      Norm1(inoi,i,1) = sum( a.hist_a(inoi,i,:) );
      if Norm1(inoi,i,1) ~= 0
         a.hist_a(inoi,i,:) = a.hist_a(inoi,i,:) / sum( a.hist_a(inoi,i,:) );
         Norm1(inoi,i,1) = sum( a.hist_a(inoi,i,:) );
      end
      Norm1(inoi,i,2) = sum( b.hist_a(inoi,i,:) );
      if Norm1(inoi,i,2) ~= 0
         b.hist_a(inoi,i,:) = b.hist_a(inoi,i,:) / sum( b.hist_a(inoi,i,:) );
         Norm1(inoi,i,2) = sum( b.hist_a(inoi,i,:) );
      end
      
      Norm1(inoi,i,3) = sum( a.hist_b(inoi,i,:) );
      if Norm1(inoi,i,3) ~= 0
         a.hist_b(inoi,i,:) = a.hist_b(inoi,i,:) / sum( a.hist_b(inoi,i,:) );
         Norm1(inoi,i,3) = sum( a.hist_b(inoi,i,:) );
      end
      Norm1(inoi,i,4) = sum( b.hist_b(inoi,i,:) );
      if Norm1(inoi,i,4) ~= 0
         b.hist_b(inoi,i,:) = b.hist_b(inoi,i,:) / sum( b.hist_b(inoi,i,:) );
         Norm1(inoi,i,4) = sum( b.hist_b(inoi,i,:) );
      end
      
      Norm1(inoi,i,5) = sum( a.hist_ax(inoi,i,:) );
      if Norm1(inoi,i,5) ~= 0
         a.hist_ax(inoi,i,:) = a.hist_ax(inoi,i,:) / sum( a.hist_ax(inoi,i,:) );
         Norm1(inoi,i,5) = sum( a.hist_ax(inoi,i,:) );
      end
      Norm1(inoi,i,6) = sum( b.hist_ax(inoi,i,:) );
      if Norm1(inoi,i,6) ~= 0
         b.hist_ax(inoi,i,:) = b.hist_ax(inoi,i,:) / sum( b.hist_ax(inoi,i,:) );
         Norm1(inoi,i,6) = sum( b.hist_ax(inoi,i,:) );
      end
      
      if i<N
         Norm1(inoi,i,7) = sum( a.hist_tri(inoi,i,:) );
         if Norm1(inoi,i,7) ~= 0
            a.hist_tri(inoi,i,:) = a.hist_tri(inoi,i,:) / sum( a.hist_tri(inoi,i,:) );
            Norm1(inoi,i,7) = sum( a.hist_tri(inoi,i,:) );
         end
         Norm1(inoi,i,8) = sum( b.hist_tri(inoi,i,:) );
         if Norm1(inoi,i,8) ~= 0
            b.hist_tri(inoi,i,:) = b.hist_tri(inoi,i,:) / sum( b.hist_tri(inoi,i,:) );
            Norm1(inoi,i,8) = sum( b.hist_tri(inoi,i,:) );
         end
      end
      
   end
end

% 2) Поиск порога
for inoi = 1:noiN
   for i = 1:N
      [~,c.T_a(inoi,i)] = max( a.hist_a(inoi,i,:) + b.hist_a(inoi,i,:) );
      [~,c.T_b(inoi,i)] = max( a.hist_b(inoi,i,:) + b.hist_b(inoi,i,:) );
      [~,c.T_ax(inoi,i)] = max( a.hist_ax(inoi,i,:) + b.hist_ax(inoi,i,:) );
      if i<N
         [~,c.T_tri(inoi,i)] = max( a.hist_tri(inoi,i,:) + b.hist_tri(inoi,i,:) );
      end
   end
end

% 3) Взаимные вероятности
for inoi = 1:noiN
   for i = 1:N
      for j = 1:ordN
         
         c.P1_a(inoi,i,j) = a.hist_a(inoi,i,j)./(a.hist_a(inoi,i,j) + b.hist_a(inoi,i,j));
         if isnan(c.P1_a(inoi,i,j))
            if j < c.T_a(inoi,i)
               c.P1_a(inoi,i,j) = 0;
            else
               c.P1_a(inoi,i,j) = 1;
            end
         end
         c.P2_a(inoi,i,j) = b.hist_a(inoi,i,j)./(a.hist_a(inoi,i,j) + b.hist_a(inoi,i,j));
         if isnan(c.P2_a(inoi,i,j))
            if j < c.T_a(inoi,i)
               c.P2_a(inoi,i,j) = 1;
            else
               c.P2_a(inoi,i,j) = 0;
            end
         end

         c.P1_b(inoi,i,j) = a.hist_b(inoi,i,j)./(a.hist_b(inoi,i,j) + b.hist_b(inoi,i,j));
         if isnan(c.P1_b(inoi,i,j))
            if j < c.T_b(inoi,i)
               c.P1_b(inoi,i,j) = 1;
            else
               c.P1_b(inoi,i,j) = 0;
            end
         end
         c.P2_b(inoi,i,j) = b.hist_b(inoi,i,j)./(b.hist_b(inoi,i,j) + a.hist_b(inoi,i,j));
         if isnan(c.P2_b(inoi,i,j))
            if j < c.T_b(inoi,i)
               c.P2_b(inoi,i,j) = 1;
            else
               c.P2_b(inoi,i,j) = 0;
            end
         end

         c.P1_ax(inoi,i,j) = a.hist_ax(inoi,i,j)./(a.hist_ax(inoi,i,j) + b.hist_ax(inoi,i,j));
         if isnan(c.P1_ax(inoi,i,j))
            if j < c.T_ax(inoi,i)
               c.P1_ax(inoi,i,j) = 0;
            else
               c.P1_ax(inoi,i,j) = 1;
            end
         end
         c.P2_ax(inoi,i,j) = b.hist_ax(inoi,i,j)./(a.hist_ax(inoi,i,j) + b.hist_ax(inoi,i,j));
         if isnan(c.P2_ax(inoi,i,j))
            if j < c.T_ax(inoi,i)
               c.P2_ax(inoi,i,j) = 1;
            else
               c.P2_ax(inoi,i,j) = 0;
            end
         end

         if i==N,break;end
         c.P1_tri(inoi,i,j) = a.hist_tri(inoi,i,j)./(a.hist_tri(inoi,i,j) + b.hist_tri(inoi,i,j));
         if isnan(c.P1_tri(inoi,i,j))
            if j < c.T_tri(inoi,i)
               c.P1_tri(inoi,i,j) = 0;
            else
               c.P1_tri(inoi,i,j) = 1;
            end
         end
         c.P2_tri(inoi,i,j) = b.hist_tri(inoi,i,j)./(a.hist_tri(inoi,i,j) + b.hist_tri(inoi,i,j));
         if isnan(c.P2_tri(inoi,i,j))
            if j < c.T_tri(inoi,i)
               c.P2_tri(inoi,i,j) = 0;
            else
               c.P2_tri(inoi,i,j) = 1;
            end
         end
         
      end
   end
end

% 4) Вероятности ошибок
for inoi = 1:noiN
   for i = 1:N
      
      c.p11_a(inoi,i) = sum( c.P1_a( inoi, i, c.T_a(inoi,i):end ) );
      c.p00_a(inoi,i) = sum( c.P2_a( inoi, i, 1:c.T_a(inoi,i) ) );
      c.p10_a(inoi,i) = sum( c.P1_a( inoi, i, 1:c.T_a(inoi,i) ) );
      c.p01_a(inoi,i) = sum( c.P2_a( inoi, i, c.T_a(inoi,i):end ) );

      c.p11_b(inoi,i) = sum( c.P1_b( inoi, i, c.T_b(inoi,i):end ) );
      c.p00_b(inoi,i) = sum( c.P2_b( inoi, i, 1:c.T_b(inoi,i) ) );
      c.p10_b(inoi,i) = sum( c.P1_b( inoi, i, 1:c.T_b(inoi,i) ) );
      c.p01_b(inoi,i) = sum( c.P2_b( inoi, i, c.T_b(inoi,i):end ) );

      c.p11_ax(inoi,i) = sum( c.P1_ax( inoi, i, c.T_ax(inoi,i):end ) );
      c.p00_ax(inoi,i) = sum( c.P2_ax( inoi, i, 1:c.T_ax(inoi,i) ) );
      c.p10_ax(inoi,i) = sum( c.P1_ax( inoi, i, 1:c.T_ax(inoi,i) ) );
      c.p01_ax(inoi,i) = sum( c.P2_ax( inoi, i, c.T_ax(inoi,i):end ) );
      
      if i<N
         c.p11_tri(inoi,i) = sum( c.P1_tri( inoi, i, c.T_tri(inoi,i):end ) );
         c.p00_tri(inoi,i) = sum( c.P2_tri( inoi, i, 1:c.T_tri(inoi,i) ) );
         c.p10_tri(inoi,i) = sum( c.P1_tri( inoi, i, 1:c.T_tri(inoi,i) ) );
         c.p01_tri(inoi,i) = sum( c.P2_tri( inoi, i, c.T_tri(inoi,i):end ) );
      end
      
   end
end

% 5) Нормировка значений ( 00+01=1  10+11=1 )
for inoi = 1:noiN
   for i = 1:N
      s = sum( c.p00_a(inoi,i) + c.p10_a(inoi,i) );   % Incorrect!!!
      c.p00_a(inoi,i) = c.p00_a(inoi,i) / s;
      c.p10_a(inoi,i) = c.p10_a(inoi,i) / s;
      s = sum( c.p01_a(inoi,i) + c.p11_a(inoi,i) );
      c.p01_a(inoi,i) = c.p01_a(inoi,i) / s;
      c.p11_a(inoi,i) = c.p11_a(inoi,i) / s;

      s = sum( c.p00_b(inoi,i) + c.p10_b(inoi,i) );
      c.p00_b(inoi,i) = c.p00_b(inoi,i) / s;
      c.p10_b(inoi,i) = c.p10_b(inoi,i) / s;
      s = sum( c.p01_b(inoi,i) + c.p11_b(inoi,i) );
      c.p01_b(inoi,i) = c.p01_b(inoi,i) / s;
      c.p11_b(inoi,i) = c.p11_b(inoi,i) / s;
      
      s = sum( c.p00_ax(inoi,i) + c.p10_ax(inoi,i) );
      c.p00_ax(inoi,i) = c.p00_ax(inoi,i) / s;
      c.p10_ax(inoi,i) = c.p10_ax(inoi,i) / s;
      s = sum( c.p01_ax(inoi,i) + c.p11_ax(inoi,i) );
      c.p01_ax(inoi,i) = c.p01_ax(inoi,i) / s;
      c.p11_ax(inoi,i) = c.p11_ax(inoi,i) / s;

      if i<N
         s = sum( c.p00_tri(inoi,i) + c.p10_tri(inoi,i) );
         c.p00_tri(inoi,i) = c.p00_tri(inoi,i) / s;
         c.p10_tri(inoi,i) = c.p10_tri(inoi,i) / s;
         s = sum( c.p01_tri(inoi,i) + c.p11_tri(inoi,i) );
         c.p01_tri(inoi,i) = c.p01_tri(inoi,i) / s;
         c.p11_tri(inoi,i) = c.p11_tri(inoi,i) / s;
      end
   end
end











