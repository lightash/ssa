function out = transform(in, to_what)
%TRANSFORM   Transforms or copies vectors or matrices between
%   [ out( 1|N , N ) ] = imp_OSR( in( 1|N , N ) , 'to_what' )
%   'to_what' can be:
%   'vector' (by default)
%   'matrix'
%   'vector_repeat'
%   'matrix_repeat'

if nargin < 2, to_what = 'vector'; end  % by default or without 2nd arg

if strcmp(to_what,'vector')
    out = in(1,:);
    for i = 2:length(in)
        out = [out in(i,:)];
    end

elseif strcmp(to_what,'matrix')
    N = floor(sqrt(length(in)));
    for i = 1:N
        out(i,:) = in( (i-1)*N+1 : i*N );
    end

elseif strcmp(to_what,'vector_repeat')
    out = in;
    for i = 2:length(in)
        out = [out in];
    end

elseif strcmp(to_what,'matrix_repeat')
    out = in;
    for i = 2:length(in)
        out = [out ; in];
    end

end

