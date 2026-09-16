%% Gauss-Hermite nodes and weights
% ref: http://www.aims.ac.za/resources/courses/na/gauss_quad.pdf
% n: order of the hermite polynomila
% x : nodes
% w : weights
function [n,x,w] = gqzero(n)
  d = sqrt((1:n)/2);
  T = diag(d,1)+diag(d,-1);
  [V,D] = eig(T);
  w = V(1,:).^2;
  w = w*sqrt(pi);
  x = diag(D);
  [x,ind] = sort(x);
  w = w(ind);
  w=w';
end

