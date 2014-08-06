clc;
close all;
clear all;

% Картинка
img = imread('Images\img0.png');
N = 350;

color = zeros(N,N,3);
for c = 1:3
    
    colstr = ['red  ';'green';'blue '];
    disp(['Processing ' colstr(c,:)])
    
    grayimg = img(:,:,c);
    % grayimg = rgb2gray(img);
    n = 1;
    a = double(grayimg(1:n:n*N,1:n:n*N));
    data = a;

    AX = floor(N/16);  % Number of axes to use

    % Sequence to choose axes
    seq = round(linspace(1,N,AX));  % Odd ones
    % seq = 1:AX;  % First ones

    % Compressing
    axis = zeros(AX,N);
    proj = axis;
    vproj = zeros(N);
    for axi = 1:AX

        axis(axi,:) = data(seq(axi),:);
        axis(axi,:) = axis(axi,:) / sqrt( axis(axi,:)*axis(axi,:)' );

        for i = 1:N
            proj(axi,i) = data(i,:) * axis(axi,:)';  % Lenghts of projections on chosen axes
            vproj(i,:) = proj(axi,i) * axis(axi,:);  % Vectors of projections
        end

        data = data - vproj;

    end

    % Decompressing
    res = zeros(N);
    for i = 1:N
        for axi = 1:AX
            res(i,:) = res(i,:) + proj(axi,i)*axis(axi,:);
        end
    end

    color(:,:,c) = res;

end

figure
subplot(131),imshow(uint8(img(1:N,1:N,:))),grid on
    title(...
    ['source (' num2str(size(img,1)) 'x' num2str(size(img,2)) 'x' num2str(size(img,3)) ')'])
subplot(132),imshow(uint8(color)),grid on
    title({
        ['result from ' num2str(AX) ' axes']
        [num2str( size(a,1)/(size(axis,1)+size(proj,1)) ,3) ' times compressed'] })
subplot(133),imshow(img(1:N,1:N,:) - uint8(color)),title('difference'),grid on


