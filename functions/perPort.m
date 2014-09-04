function [portrait] = perPort(input,T)
%PERPORT   Periodic portrait of a signal
%   input - signal vector
%   T - vector of periods to compute

% through periods
for per = 1:length(T)
   N = T(per);
   N1 = N^2;
   
   % through windows
   for win = 1:length(input) - (N1-1)
      
      data_win = input( win : win + (N1-1) );
      
      % pre-offsetting data for phase correction
      ind = round( win-1 - N1*floor( (win-1)/N1 ) );  % globally (from end to start)
      data_offset = [data_win(N1+1-ind:N1) data_win(1:N1-ind)];
      
      % through N in ring (window end is connected to start)
      data_win_mat = transform(data_offset,'matrix');
      [vproj, portrait.period{per}.window{win}.win_basis] = GSOrth(imp_OSR(data_win_mat));
      portrait.period{per}.window{win}.svproj(1,:) = sum(vproj,2);
      
      for rin = 2:N
         
         data_ring = transform( [data_win(N1+1-rin:N1) data_win(1:N1-rin)] ,'matrix');
         
         portrait.period{per}.window{win}.svproj(rin,:) = sum( portrait.period{per}.window{win}.win_basis * data_ring' ,1)';
         
      end
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
%          portrait.period{per}.window{win}.per_offset(imp,:) = ...
%             [portrait.period{per}.window{win}.per_win( imp, N+1-ind : N ),...
%             portrait.period{per}.window{win}.per_win( imp, 1 : N-ind )];
