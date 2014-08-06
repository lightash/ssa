function [a, b, e] = ABAxTriDecomp(A, noi, stn, ord, fname)
%ABAxTriDecomp   Decomposition of matrix A to vector B and
% matrices Ax and Tri and their correlations
%   noi - noise vector
%   stn - statistics size
%   ord - positions of histogram bars
%   fname - filename to save variables

a = struct;
b = struct;
e = struct;

nosave = 0; if nargin < 5, nosave = 1; end  % by default or without 5th arg

N = size(A,1);

% Standard
[e.a2, e.b2] = OSR(A);

[~, e.ax, e.lax, e.tri] = GSOrth(e.a2);  % Orthogonalization

for axi = 1:N  % Normalization
   e.la(axi) = sqrt( A(axi,:)*A(axi,:)' );
   e.a(axi,:) = A(axi,:)/e.la(axi);
   if axi<N
      e.ltri(axi) = sqrt(e.tri(axi,:)*e.tri(axi,:)');
      e.tri(axi,:) = e.tri(axi,:)/e.ltri(axi);
   end
end

e.lb = sqrt(e.b2 * e.b2');
e.b = e.b2 / e.lb;

for i = 1:N
   e.cor_b(i) = e.a(i,:) * e.b';
end

% % ...for double components
%    for axi = 1:N-1
%       e2.la(axi) = sqrt( [A(axi,:) A(axi+1,:)] * [A(axi,:) A(axi+1,:)]' );
%       e2.a(axi,:) = [A(axi,:) A(axi+1,:)] / e2.la(axi);
%       if axi<N-1
%          e2.ltri(axi) = sqrt( [e.tri(axi,:) e.tri(axi+1,:)] * [e.tri(axi,:) e.tri(axi+1,:)]' );
%          e2.tri(axi,:) = [e.tri(axi,:) e.tri(axi+1,:)] / e2.ltri(axi);
%       end
%    end
% 
%    e2.lb = sqrt( [e.b2 e.b2] * [e.b2 e.b2]' );
%    e2.b = [e.b2 e.b2] / e2.lb;
% 
%    for i = 1:N-1
%       e2.cor_b(i) = e2.a(i,:) * e2.b';
% 
%       e2.lax(i,:) = sqrt( [e.ax(i,:) e.ax(i+1,:)] * [e.ax(i,:) e.ax(i+1,:)]' );
%       e2.ax(i,:) = [e.ax(i,:) e.ax(i+1,:)] / e2.lax(i,:);
%    end

if nargin > 1  % if correlations are needed

   % With noise
   noiN = length(noi);
   ordN = length(ord);

   disp(['Noises from ' num2str(noi(1),'%1.2f ') ' to ' num2str(noi(end),'%1.2f ') ' at ' num2str(noi(2)-noi(1),'%1.2f ')])
   disp(['   ' num2str(stn) ' times each'])

   a.hist_a = zeros(noiN,N,ordN);
   a.hist_b = zeros(noiN,N,ordN);
   a.hist_ax = zeros(noiN,N,ordN);
   a.hist_tri = zeros(noiN,N,ordN);

   b.hist_a = zeros(noiN,N,ordN);
   b.hist_b = zeros(noiN,N,ordN);
   b.hist_ax = zeros(noiN,N,ordN);
   b.hist_tri = zeros(noiN,N,ordN);

   % % ...for double components
   %    a2.hist_a = zeros(noiN,N,ordN);
   %    a2.hist_b = zeros(noiN,N,ordN);
   %    a2.hist_ax = zeros(noiN,N,ordN);
   %    a2.hist_tri = zeros(noiN,N,ordN);
   % 
   %    b2.hist_a = zeros(noiN,N,ordN);
   %    b2.hist_b = zeros(noiN,N,ordN);
   %    b2.hist_ax = zeros(noiN,N,ordN);
   %    b2.hist_tri = zeros(noiN,N,ordN);

   for inoi = 1:noiN   

      % Bar in console
      numBars = 50;
      thrBar = inoi/noiN * numBars;
      calc = '';
      for i = 1:numBars
         if i <= thrBar
            calc = [calc '#'];
         else
            calc = [calc '-'];
         end
      end
      disp([calc '    ' num2str( 100*inoi/noiN, '%3.1f' ) '%    Noise = ' num2str( noi(inoi) )])


      % Gathering statistics
      for stat = 1:stn

         noise = noi(inoi) * randn(N);
         a.a = A + noise;


         % For signal+noise
         [a.a2,a.b2] = OSR(a.a);

         [~, a.ax, a.lax, a.tri] = GSOrth(a.a2);  % Orthogonalization

         for i = 1:N  % Normalization
            if i<N
               a.ltri(i) = sqrt(a.tri(i,:)*a.tri(i,:)');
               a.tri(i,:) = a.tri(i,:)/a.ltri(i);
            end
         end

   %       % ...for double components
   %          for axi = 1:N-2
   %             a2.ltri(axi) = sqrt( [a.tri(axi,:) a.tri(axi+1,:)] * [a.tri(axi,:) a.tri(axi+1,:)]' );
   %             a2.tri(axi,:) = [a.tri(axi,:) a.tri(axi+1,:)] / a2.ltri(axi);
   %          end

         for i = 1:N  % Correlations

            a.cor_a(stat,i) = e.a(i,:) * ( a.a(i,:)/sqrt( a.a(i,:)*a.a(i,:)' ) )';
            a.cor_b(stat,i) = e.a(i,:) * ( a.b2/sqrt( a.b2*a.b2' ) )';
            a.cor_ax(stat,i) = e.ax(i,:) * a.ax(i,:)';

            if i==N,break;end
            a.cor_tri(stat,i) = e.tri(i,:) * a.tri(i,:)';
         end

   %       % ...for double components
   %          for i = 1:N-1
   %             a2.cor_a(stat,i) = e2.a(i,:) * ( [a.a(i,:) a.a(i+1,:)] / sqrt( [a.a(i,:) a.a(i+1,:)] * [a.a(i,:) a.a(i+1,:)]' ) )';
   %             a2.cor_b(stat,i) = e2.a(i,:) * ( [a.b2 a.b2] / sqrt( [a.b2 a.b2] * [a.b2 a.b2]' ) )';
   %             a2.cor_ax(stat,i) = e2.ax(i,:) * [a.ax(i,:) a.ax(i+1,:)]';
   %             if i<N-1
   %                a2.cor_tri(stat,i) = e2.tri(i,:) * [a.tri(i,:) a.tri(i+1,:)]';
   %             end
   %          end


         % For noise
         [b.a2,b.b2] = OSR(noise);

         [~, b.ax, b.lax, b.tri] = GSOrth(b.a2);  % Orthogonalization

         for axi = 1:N  % Normalization
            if axi<N
               b.ltri(axi) = sqrt(b.tri(axi,:)*b.tri(axi,:)');
               b.tri(axi,:) = b.tri(axi,:)/b.ltri(axi);
            end
         end

         for i = 1:N  % Correlations

            b.cor_a(stat,i) = e.a(i,:) * ( noise(i,:)/sqrt( noise(i,:)*noise(i,:)' ) )';
            b.cor_b(stat,i) = e.a(i,:) * ( b.b2/sqrt( b.b2*b.b2' ) )';
            b.cor_ax(stat,i) = e.ax(i,:) * b.ax(i,:)';

            if i==N,break;end
            b.cor_tri(stat,i) = e.tri(i,:) * b.tri(i,:)';
         end

   %       % ...for double components
   %          for i = 1:N-1
   %             b2.cor_a(stat,i) = e2.a(i,:) * ( [noise(i,:) noise(i+1,:)] / sqrt( [noise(i,:) noise(i+1,:)] * [noise(i,:) noise(i+1,:)]' ) )';
   %             b2.cor_b(stat,i) = e2.a(i,:) * ( [b.b2 b.b2] / sqrt( [b.b2 b.b2] * [b.b2 b.b2]' ) )';
   %             b2.cor_ax(stat,i) = e2.ax(i,:) * [b.ax(i,:) b.ax(i+1,:)]';
   %             if i<N-1
   %                b2.cor_tri(stat,i) = e2.tri(i,:) * [b.tri(i,:) b.tri(i+1,:)]';
   %             end
   %          end

      end

      for i = 1:N  % Histograms for statistics st_hist -> hist

         a.hist_a(inoi,i,:) = hist(a.cor_a(:,i),ord);
         b.hist_a(inoi,i,:) = hist(b.cor_a(:,i),ord);

         a.hist_b(inoi,i,:) = hist(a.cor_b(:,i),ord);
         b.hist_b(inoi,i,:) = hist(b.cor_b(:,i),ord);

         a.hist_ax(inoi,i,:) = hist(a.cor_ax(:,i),ord);
         b.hist_ax(inoi,i,:) = hist(b.cor_ax(:,i),ord);

         if i==N,break;end
         a.hist_tri(inoi,i,:) = hist(a.cor_tri(:,i),ord);
         b.hist_tri(inoi,i,:) = hist(b.cor_tri(:,i),ord);
      end

   %    % ...for double components
   %       for i = 1:N-1
   %          a2.hist_a(inoi,i,:) = hist(a2.cor_a(:,i),ord);
   %          b2.hist_a(inoi,i,:) = hist(b2.cor_a(:,i),ord);
   %          
   %          a2.hist_b(inoi,i,:) = hist(a2.cor_b(:,i),ord);
   %          b2.hist_b(inoi,i,:) = hist(b2.cor_b(:,i),ord);
   % 
   %          a2.hist_ax(inoi,i,:) = hist(a2.cor_ax(:,i),ord);
   %          b2.hist_ax(inoi,i,:) = hist(b2.cor_ax(:,i),ord);
   % 
   %          if i==N-1,break;end
   %          a2.hist_tri(inoi,i,:) = hist(a2.cor_tri(:,i),ord);
   %          b2.hist_tri(inoi,i,:) = hist(b2.cor_tri(:,i),ord);
   %       end


      a.am_cor_a(inoi,:) = abs(mean(a.cor_a,1));
      a.am_cor_b(inoi,:) = abs(mean(a.cor_b,1));
      a.am_cor_ax(inoi,:) = abs(mean(a.cor_ax,1));
      a.am_cor_tri(inoi,:) = abs(mean(a.cor_tri,1));
      b.am_cor_a(inoi,:) = abs(mean(b.cor_a,1));
      b.am_cor_b(inoi,:) = abs(mean(b.cor_b,1));
      b.am_cor_ax(inoi,:) = abs(mean(b.cor_ax,1));
      b.am_cor_tri(inoi,:) = abs(mean(b.cor_tri,1));

   %    % ...for double components
   %       a2.am_cor_a(inoi,:) = abs(mean(a2.cor_a,1));
   %       a2.am_cor_b(inoi,:) = abs(mean(a2.cor_b,1));
   %       a2.am_cor_ax(inoi,:) = abs(mean(a2.cor_ax,1));
   %       a2.am_cor_tri(inoi,:) = abs(mean(a2.cor_tri,1));
   %       b2.am_cor_a(inoi,:) = abs(mean(b2.cor_a,1));
   %       b2.am_cor_b(inoi,:) = abs(mean(b2.cor_b,1));
   %       b2.am_cor_ax(inoi,:) = abs(mean(b2.cor_ax,1));
   %       b2.am_cor_tri(inoi,:) = abs(mean(b2.cor_tri,1));

   end

   if nosave == 0
      save(fname);
   end

end