function res = u0(x, v1, n, p1)
% gives the initial profile: first n elements are the first derivative,
% then the ball's first derivative. next n elements are the positions, then
% the ball's position

x = [x;0;0;0];
res0 = 0.*repmat(x,2,1);
res0(n+4) = v1;
res0(end) = p1;

res = res0;
end
