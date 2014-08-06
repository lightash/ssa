function [B, A3] = imp_OSR(a1,path)
%IMP_OSR   Continuous and impulse components of a square matrix
%   [ B(N_imp,N), A3(N,N) ] = imp_OSR( A1(N,N) ) returns the 
%   continuous and impulse components of a square matrix A1
%   and remains from crossings out.

if nargin < 2, path = 'from_end'; end  % by default or without 2nd arg

N = length(a1);

A1 = zeros(N,N,N);
A2 = zeros(N,N,N);
A3 = zeros(N,N,N);
B = zeros(N);

A1(:,:,1) = a1;

[B(1,:), A2(:,:,1)] = OSR(A1(:,:,1));  % Continuous component

IJ = 1:N;  % List of used rows & columns

for imp = 2:N  % Impulse component
    
    if strcmp(path,'energy')  % Path of max energy
        imp_en = zeros(1,N+2-imp);
        for rc = 1 : N+2-imp  % Rows & columns
            choise_ind = IJ(1:rc-1);
            choise_ind = [choise_ind IJ(rc+1:end)]; %#ok<AGROW>

            b = OSR( A2(choise_ind, choise_ind, imp) );
            imp_en(rc) = b*b';  % Energies of components
        end
        [~,next] = max(imp_en);
    end
    
    if strcmp(path,'from_end')  % Crossing out from end
        next = N-imp+2;
    end

    IJ = [IJ(1:next-1) IJ(next+1:end)];
    
    A3(next,IJ,imp) = A2(next,IJ,imp-1);
    A3(IJ,next,imp) = A2(IJ,next,imp-1);
    A3(next,next,imp) = A2(next,next,imp-1);
    
    A1(IJ,IJ,imp) = A2(IJ,IJ,imp-1);
    
    [B(imp,IJ), A2(IJ,IJ,imp)] = OSR(A1(IJ,IJ,imp));
    
end