function [vproj, axes, lax, tri] = GSOrth(in)
%GSORTH   Gram-Shmidt orthogonalization
%   vproj - vectors of projections on chosen axes
%   axes - chosen axes
%   lax - lengths of vproj
%   tri - remainings from orthogonalizing,
%      triangular matrix

N = size(in,1);

vproj = zeros(N);
axes = zeros(N);
lax = zeros(1,N);
tri = zeros(N);

for axi = 1:N

   % Normalization
   lax(axi) = len(in(axi,:));
   axes(axi,:) = in(axi,:) / lax(axi);

   for i = axi+1:N
      tri(axi,i) = axes(axi,:) * in(i,:)';  % Lenghts of projections on chosen axes
      vproj(i,:) = tri(axi,i) * axes(axi,:);  % Vectors of projections
      in(i,:) = in(i,:) - vproj(i,:);
   end
end

vproj = in;