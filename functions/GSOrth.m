function [vproj, axes, lax, other] = GSOrth(data)
%GSORTH   Gram-Shmidt orthogonalization
%   vproj - vectors of projections on chosen axes
%   axes - chosen axes
%   lax - lengths of vproj
%   other - remainings from orthogonalizing,
%      triangular matrix

N = size(data,1);

vproj = zeros(N);
axes = zeros(N);
lax = zeros(1,N);
other = zeros(N);
proj = zeros(N);

for axi = 1:N

   % Normalization
   axes(axi,:) = data(axi,:);
   lax(axi) = sqrt( axes(axi,:)*axes(axi,:)' );
   axes(axi,:) = axes(axi,:) / lax(axi);

   for i = axi+1:N
      proj(axi,i) = axes(axi,:) * data(i,:)';  % Lenghts of projections on chosen axes
      vproj(i,:) = proj(axi,i) * axes(axi,:);  % Vectors of projections
      data(i,:) = data(i,:) - vproj(i,:);
      other(axi,i) = proj(axi,i);
   end
end

vproj = data;