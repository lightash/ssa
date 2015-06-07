function [E, C, A3] = impAM(a1, path)
%impAM   Continuous and impulse components of a matrix.
%   [ E(N_imp,N), C(N_imp,M), A3(M,N) ] = impAM( A1(M,N) ) returns the 
%   continuous ( E(1,:) ) and impulse ( E(2:end, :) ) components and
%   their amplitudes along axes ( C ) of a matrix A1 
%   and remains from crossings out ( A3 ).
%   path - string that defines rows & columns that are crossed out:
%   'from_end' - crossing out from end;
%   'energy' - path of max energy.

if nargin < 2, path = 'energy'; end  % by default

M = size(a1,1);
N = M;% TODO: size(a1,2);

A1 = zeros(M,N,M);
A2 = zeros(M,N,M);
A3 = zeros(M,N);
E = zeros(M,N);
C = zeros(N,M);

A1(:,:,1) = a1;

[E(1,:), C(1,:), A2(:,:,1)] = AM( A1(:,:,1) );  % Continuous component

UR = 1:M;  % List of used rows
UC = 1:N;  % List of used columns

for imp = 2:N  % Impulse component
   
   % Path of max energy
   if strcmp(path,'energy')

      imp_len = zeros(M - imp+2, N - imp+2);
      for row = 1 : M - imp+2
         choice_row = UR(1:row-1);
         choice_row = [choice_row UR(row+1:end)]; %#ok<AGROW>

         for col = row% TODO: 1:N - imp+2
            choice_col = UC(1:col-1);
            choice_col = [choice_col UC(col+1:end)]; %#ok<AGROW>

            [e, c] = AM( A2(choice_row, choice_col, imp-1) );
            e = e*mean(c);
            imp_len(row, col) = e*e';  % Energies of components
         end
      end

      [max_row, next_row] = max(imp_len);
      [~,next(2)] = max(max_row);
      next(1) = next_row(next(2));
   end

   % Crossing out from end
   if strcmp(path,'from_end')
      next = [M - imp+2, N - imp+2];
   end

   UR = [UR(1:next(1)-1) UR(next(1)+1:end)];
   UC = [UC(1:next(2)-1) UC(next(2)+1:end)];

   A3(next(1),UC) = A3(next(1),UC) + A2(next(1),UC,imp-1);
   A3(UR,next(2)) = A3(UR,next(2)) + A2(UR,next(2),imp-1);
   A3(next(1),next(2)) = A3(next(1),next(2)) + A2(next(1),next(2),imp-1);

   A1(UR,UC,imp) = A2(UR,UC,imp-1);

   [E(imp,UC), C(imp,UR), A2(UR,UC,imp)] = AM( A1(UR,UC,imp) );

end