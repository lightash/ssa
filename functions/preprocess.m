function [sr_mat, a_coll, a_ort, w] = preprocess(a)
%PREPROCESS   Implements the modified version of OSR
%   [sr_mat(N,N), a_coll(N,N), a_ort(N,N), w(N)] = preprocess(a(N,N))
%   returns the mean-in-matrix component of a square matrix A1
%   along with its collinear a_coll, orthogonal a_ort components
%   and lengthes of collinear vectors w.
    
    N = size(a,1);
    
    sr = mean(a);               % Среднее значение строк
    sr_mat = transform(sr,'matrix_repeat');
    
    len_e_sr = sqrt(sr*sr');    % Длина орта среднего
    e_sr = sr./len_e_sr;        % Орт среднего

    a_ost = zeros(N);           % Выделение памяти
    w = zeros(1,N);
    a_coll = zeros(N);
    a_ort = zeros(N);
    
    for i = 1:N
        a_ost(i,:) = a(i,:) - sr;               % Остатки строк от среднего
        w(i) = e_sr * a_ost(i,:)';              % Проекции остатков на орт
        a_coll(i,:) = w(i) * e_sr;              % Коллинеарные орту остатки
        a_ort(i,:) = a_ost(i,:) - a_coll(i,:);  % Ортогональные орту остатки
    end
    
end