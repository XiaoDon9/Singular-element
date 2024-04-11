clear
clc

syms r s

% h(8) = 1/2 * (1 - s^2) * (1 + r);
% h(7) = 1/2 * (1 - r^2) * (1 - s);
% h(6) = 1/2 * (1 - s^2) * (1 - r);
% h(5) = 1/2 * (1 - r^2) * (1 + s);
%
% h(4) = 1/4 * (1 + r) * (1 - s) - 1/2 * h(7) - 1/2 * h(8);
% h(3) = 1/4 * (1 - r) * (1 - s) - 1/2 * h(6) - 1/2 * h(7);
% h(2) = 1/4 * (1 - r) * (1 + s) - 1/2 * h(5) - 1/2 * h(6);
% h(1) = 1/4 * (1 + r) * (1 + s) - 1/2 * h(5) - 1/2 * h(8);


% simplified_h = simplify(h);

cord_loc = [[1, -1, -1, 1, 0, -1, 0, 1].', [1, 1, -1, -1, 1, 0, -1, 0].'];

for i = 1:8
    h(i) = Ni(r, s, cord_loc(i, :));
end

simplified_h = simplify(h);

h236 = h(2)+ h(3) +h(6);
simplified_h236 = simplify(h236);

% simplify(h - h2)
% % ans
% % [0, 0, 0, 0, 0, 0, 0, 0]

cord_glb = [[2, 0, 0, 2, 0.5, 0, 0.5, 2].', [1, 0, 0, -1, 0.25, 0, -0.25, 0].'];

X = cord_glb([1,4,5,7,8],:)' * simplified_h([1,4,5,7,8]).' +cord_glb(2,:)'* simplified_h236 ;
simplify_X = simplify(X);

substitutedExpr1 = simplify(subs(X, {s}, {0})); % 自然坐标与总体坐标的关系 ！看阶次
% r+1 = (2*x)^(1/2)
clear X i h

J = [diff(simplify_X, r), diff(simplify_X, s)];
J_inv = inv(J);

substitutedExpr2 = simplify(subs(J_inv, {s}, {0})); % 雅可比矩阵的奇异性 ！初步判断

dNdr = diff([simplified_h([1,4,5,7,8]),simplified_h236], r);
dNds = diff([simplified_h([1,4,5,7,8]),simplified_h236], s);

DNDx = J \ [dNdr; dNds];
clear dNdr dNds

substitutedExpr3 = simplify(subs(DNDx, {s}, {0})); % 应变和自然坐标之间的关系 ！看阶次
% 上式化简后，各项分母只含（1+r）的一次项
% 所以是只有R^(-1/2)的奇异阶

function Ni = Ni(r, s, cord)
    ri = cord(1);
    si = cord(2);
    Ni = ((1 + r * ri) * (1 + s * si) - (1 - r ^ 2) * (1 + s * si) - ...
        (1 - s ^ 2) * (1 + r * ri)) * ri ^ 2 * si ^ 2/4 + ...
        (1 - r ^ 2) * (1 + s * si) * (1 - ri ^ 2) * si ^ 2/2 + ...
        (1 - s ^ 2) * (1 + r * ri) * (1 - si ^ 2) * ri ^ 2/2;
end
