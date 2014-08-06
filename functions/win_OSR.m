function [W, WE, A3] = win_OSR(data, N)
%WIN_OSR   Shifted periodic components of a signal
%   [ W(Shifts+1,N), WE(1,Shifts+1), A3(N,N,Shifts+1) ] =
%   = win_OSR(data, N) returns the shifted periodical
%   components of a signal data W,
%   their energy WE and remains A3.

N1 = N^2;
Shifts = length(data) - N^2;  % Number of shifts

W = zeros(Shifts+1,N);
WE = zeros(1,Shifts+1);
A3 = zeros(N,N,Shifts+1);

for win = 1:Shifts+1  % Shifted windows
%     disp(['    window = ' num2str(win)])
    
    A1 = transform(data( win : win-1 + N1 ),'matrix');  % Shifting
    
    [B, ~] = imp_OSR(A1);  % ~ -> A3(:,:,win)
    
    window_sum = sum(B,2)';
    
    if win > N
        ind = win - N*floor((win-1)/N);
    else
        ind = win;
    end
    
%     W(win,:) = window_sum;
    W(win,:) = [window_sum(N+2-ind:N) window_sum(1:N+1-ind)];  % Windows with shifts
    WE(win) = window_sum * window_sum';  % Energies of windows
end