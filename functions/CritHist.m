function [H] = CritHist(S)
% Histograms for formal criteriums.
%    S is a structure of samples of copeting standards;
%    H is a structure with histograms of samples.

winL = size(S{1},2);
Nbins = 16;
sampL = [size(S{1},1) size(S{2},1)];

scale = zeros(winL, Nbins);
for i = 1:winL
   [~,scale(i,:)] = hist( [S{1};S{2}], Nbins);
   
   for type = 1:2
      H{type}(:,i) = hist( S{type}(:,i) ,scale(i,:)) / sampL(type)';
      % Распределить наименьший столбец по нулям
      minh = min(H{type}( H{type}(:,i)>0, i));
      mzh = minh/( length( H{type}( H{type}(:,i)==0, i)) +1);
      H{type}( H{type}(:,i)==0 | H{type}(:,i)==minh, i) = mzh;
   end
end
