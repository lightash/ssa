function [b2, a2,D,Q] = OSR(a1)
%OSR   Periodic component of a square matrix
%   [ B2(N), A2(N,N) ] = OSR( A1(N,N) ) returns the 
%   periodic component of a square matrix A1
%   and remains from it.

% if nargin < 2, out = 0; end

N = size(a1,1);

%------------------------------------#1------------------------------------
% ������������ ��������� �������������� ��� 1 ������

c = zeros(N);

for i = 2:N                        % ������
   
   c(i,:) = a1(i,:) - a1(1,:);
   
end

% %-----------------------------------#1.5-----------------------------------
% % ���������� ����� c[N][N]
% 
% for i = 2:N
%    c(i,:) = c(i,:) / max(abs(c(i,:)));
% end

%------------------------------------#2------------------------------------
% ������ D[j] - ��������������� ���������� j-�� �������� 1-�� ������

D = zeros(1,N);

% cc = zeros(N-1);
% f = 1;
% for j = 1:N                    % �� ��������� ����������
%     for i = 1:N-1              % �� ������� �� ������
%         for k = 1:N-1          % �� �������� ����� �������
%             if k >= j
%                 cc(i,k) = c(i+1,k+1);
%             else
%                 cc(i,k) = c(i+1,k);
%             end
%         end
%     end
%     D(j) = det(cc) * f;
%     f = -f;
% end

for j = 1:N  % �� ��������� ����������
   cc = [c( 2:N, 1:j-1 ) c( 2:N, j+1:N )];
   D(j) = det(cc) * (-1)^(j-1);
end

%-----------------------------------#3,4-----------------------------------
% ������ Q - ����������� �������� ����������� � ����� b2

Q = sqrt(D*D');

if Q == 0
   Q = 1e-32;
end

%------------------------------------#5------------------------------------
% ������ e2 - ��������� ����������� b2

e2 = D/Q;

%------------------------------------#6------------------------------------
% ������ p - ����� b2

p = a1(1,:)*e2';

%------------------------------------#7------------------------------------
% ������ ��������� ������� b2 - ����� ������������ ����� ������� a1

b2 = p*e2;

%------------------------------------#8------------------------------------
% ������ a2[1:N][1:N] - ������� ��������

a2 = zeros(N);

for i = 1:N
   a2(i,:) = a1(i,:) - b2;
end

% %------------------------------------#12-----------------------------------
% % ������ (b2,a2) - ���������� ������������ ��� �������� ���������������
% 
% if out == 1
%     SK = b2*a2';
%     disp(['lg(SK)_' num2str(N) ' -> ' num2str(round(log10(norm(SK,inf))))])
% end
