function [xqrt,fqrt,fq] = QRT(f,qs,smooth,integrate)
%QRT Quasi-reverse transformation
%   [ xqrt,fqrt,fq ] = QRT( x,f,qs ) returns the 
%   quasi-reverse transformation ( xqrt,fqrt ) and quanted function ( fq ) 
%   of input function ( f ) with quanting step ( qs ).
%   smooth - length of transition to be smoothed;
%   integrate - whether to inegrate instead of absolute values.

if nargin < 3, smooth = 0; end
if nargin < 4, integrate = 0; end

x = 1:length(f);
fq = qs*round( f/qs );  % Quanted function

df = [0 diff(fq)];
xt = [1 x(df ~= 0)]; % Transition grid

ft = fq(xt);         % Transition function

if smooth ~= 0
   ft = [ft(1) ft ft(end)];
   xt = [0 xt xt(end)+1];
   i = 2;
   while i<length(xt)
      if xt(i)-xt(i-1)<smooth
         if ft(i-1)==ft(i+1)
            if i>2 && i<length(xt)-1 && ft(i-2)-ft(i-1)~=ft(i+2)-ft(i+1)
               xt = [xt(1:i-1) xt(i+1:end)];
               ft = [ft(1:i-1) ft(i+1:end)];
               i = i-1;
            end
         end
      end
      i = i+1;
   end
   ft = ft(2:end-1);
   xt = xt(2:end-1);
end

xqrt = [xt length(f)+1];

if integrate
   fqrt(1) = 0;
   for i = 2:length(xt)
      fqrt(i) = fqrt(i-1) + sign(ft(i)-ft(i-1))*(xt(i)-xt(i-1));
   end
   fqrt(i+1) = length(f)-xt(end);
else
   fqrt = [0 sign(diff(ft)).*diff(xt)];  % Rotated function
   fqrt = [fqrt length(f)-xt(end)];
end

