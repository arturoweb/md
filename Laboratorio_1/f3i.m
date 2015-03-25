function [s] = f3i(x)
    s = [ x(1) - 4*x(1)^2 - x(1)*x(2),
          2*x(2) - x(2)^2 - 3*x(1)*x(2) ];
end

