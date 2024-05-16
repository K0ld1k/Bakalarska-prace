function [sigma] = sigma(z)
%  SIGMA V�po�et sigma
epsilon = 0.1;
sigma = (1/epsilon)*(sqrt(1 + (epsilon * norm(z)^2)) - 1);
end