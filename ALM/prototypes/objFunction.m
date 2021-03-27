function [f,g] = objFunction(x)

f  = (x(1) - 1)^2 + (x(2) - 2.5)^2;

if nargout > 1
  g  = [2*x(1) - 2; 2*x(2) - 5];
end
