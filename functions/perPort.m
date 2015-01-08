function [portrait] = perPort(input,T,imp_meth)
%PERPORT   Periodic portrait of a signal
%   input - signal vector
%   T - vector of periods to compute
%   imp_meth - string that defines which impulse matrix method to use

if nargin < 3, imp_meth = 'OSR'; end  % by default

% through periods
for per = 1:length(T)
   N = T(per);
   N1 = N^2;
   
   % through windows
   for win = 1:length(input) - (N1-1)
      
      data_win = input( win : win + (N1-1) );
      
      % pre-offsetting data for phase correction
%       ind = round( win-1 - N1*floor( (win-1)/N1 ) );  % globally (from end to start)
%       data_offset = [data_win(N1+1-ind:N1) data_win(1:N1-ind)];
      data_offset = data_win;
      
       data_win_mat = transform(data_offset,'matrix');
%        portrait{per}{win}.data_win_mat = data_win_mat;
      
      % Use imp_OSR function
      if strcmp(imp_meth,'OSR')
         [nonOrth, portrait{per}{win}.rem] = imp_OSR(data_win_mat,'energy');
         
      % Use impAM function
      elseif strcmp(imp_meth,'AM')
         [nonOrth, Amp, portrait{per}{win}.rem] = impAM(data_win_mat,'energy');
         MaxAmp = zeros(1,size(Amp,1));
         for i = 1:size(Amp,1)
            MaxAmp(i) = mean(Amp(i,:));
%             if MaxAmp(i) < 1e-12
%                MaxAmp(i) = min(Amp(i,:));
%             end
            nonOrth(i,:) = nonOrth(i,:) * MaxAmp(i);
         end
      end
      
      portrait{per}{win}.nonOrth = nonOrth;
      
      [vproj, portrait{per}{win}.win_basis] = GSOrth(nonOrth);
      portrait{per}{win}.svproj(1,:) = sum(vproj,1);
      
%       % through N in ring (window end is connected to start)
%       for rin = 2:N
%          
%          data_ring = transform( [data_win(N1+1-rin:N1) data_win(1:N1-rin)] ,'matrix');
%          
%          portrait{per}{win}.svproj(rin,:) = sum( portrait{per}.window{win}.win_basis * data_ring' ,1)';
%          
%       end
   end
end

%          % pre-offsetting data for phase correction
%          ind = round( win-1 - N1*floor( (win-1)/N1 ) );  % globally (from end to start)
%          data1 = [data_win(N1+1-ind:N1) data_win(1:N1-ind)];
%          data_offset = transform( data1 ,'matrix');
%          ind = round( win-1 - N*floor( (win-1)/N ) );  % locally (in each period N)
%          data1 = transform( datapart ,'matrix');
%          data_offset = [data1(:, N+1-ind:N) data1(:, 1:N-ind)];
%          
%          % post-offsetting data for same phase in each window
%          ind = round( win-1 - N*floor( (win-1)/N ) );  % locally (in each period N)
%          portrait{per}{win}.per_offset(imp,:) = ...
%             [portrait{per}{win}.per_win( imp, N+1-ind : N ),...
%             portrait{per}{win}.per_win( imp, 1 : N-ind )];
