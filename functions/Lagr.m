function [Y,A,B] = Lagr(x,y,X)

n = length(x);
N = length(X);

Y = zeros(1,N);
for i = 1:N
   for j = 1:n
      for k = 1:n
         a(k) = (X(i)-x(k))*1e-4;
         b(k) = (x(j)-x(k))*1e-4;
      end
      ac = [a(1:j-1) a(j+1:n)];
      bc = [b(1:j-1) b(j+1:n)];
      A(i,j) = prod(ac);
      B(i,j) = prod(bc);
   end
   Y(i) = sum(A(i,:)./B(i,:).*y);
end
