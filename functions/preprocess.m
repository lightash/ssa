function [sr_mat, a_coll, a_ort, w] = preprocess(a)
%PREPROCESS   Implements the modified version of OSR
%   [sr_mat(N,N), a_coll(N,N), a_ort(N,N), w(N)] = preprocess(a(N,N))
%   returns the mean-in-matrix component of a square matrix A1
%   along with its collinear a_coll, orthogonal a_ort components
%   and lengthes of collinear vectors w.
    
    N = size(a,1);
    
    sr = mean(a);               % ������� �������� �����
    sr_mat = transform(sr,'matrix_repeat');
    
    len_e_sr = sqrt(sr*sr');    % ����� ���� ��������
    e_sr = sr./len_e_sr;        % ��� ��������

    a_ost = zeros(N);           % ��������� ������
    w = zeros(1,N);
    a_coll = zeros(N);
    a_ort = zeros(N);
    
    for i = 1:N
        a_ost(i,:) = a(i,:) - sr;               % ������� ����� �� ��������
        w(i) = e_sr * a_ost(i,:)';              % �������� �������� �� ���
        a_coll(i,:) = w(i) * e_sr;              % ������������ ���� �������
        a_ort(i,:) = a_ost(i,:) - a_coll(i,:);  % ������������� ���� �������
    end
    
end