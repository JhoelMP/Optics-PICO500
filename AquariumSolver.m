function [v dvdt] = AquariumSolver(t, x, y, n)

num_layers = numel(y);
num_t = numel(t);
size_orig = size(t);

n = n(1)./n;

t = repmat(t(:),1,num_layers);
y = repmat(y(:)',num_t,1);
n = repmat(n(:)',num_t,1);


v = sum(y.*n.*t.*((1 - (n.*t).^2).^-.5), 2) - x;

v = reshape(v,size_orig);

dvdt = sum(y.*n.*( ...
    ((1 - (n.*t).^2).^-.5) + ...
    ((n.*t).^2).*((1 - (n.*t).^2).^-1.5) ), 2);

