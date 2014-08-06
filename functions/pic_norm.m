function out = pic_norm(in, type)

if nargin < 2, type = 'middle'; end  % by default or without 2nd arg

if strcmp(type,'to_one')
    
    min_in = min(min(in));
    from_zero = in - min_in;
    max_in = max(max(from_zero));

    out = from_zero / max_in;
    
elseif strcmp(type,'middle')
    
    out = in/2 + 127;
    
end