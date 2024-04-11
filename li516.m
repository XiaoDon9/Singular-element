clear
clc

syms r s t

col = [1, -1, -1, 1];
col_rep = repmat(col, 1, 2);
col2 = [1, 1, -1, -1];
col_rep2 = repmat(col2, 1, 2);
col_rep3 = [ones(1, 4), -ones(1, 4)];

cord_loc = [col_rep', col_rep2', col_rep3'];
cord_glb = cord_loc + 1;

clear col_rep3 col_rep2 col_rep col2 col

for i = 1:8
    h(i) = g(r,s,t,cord_loc(i, :));
end

h1234 = 0;

for i = 1:4
    h1234 = h1234 + h(i);
end

simplified_h1234 = simplify(h1234);
clear h1234

h56 = 0;

for i = 5:6
    h56 = h56 + h(i);
end

simplified_h56 = simplify(h56);
clear h56

X = cord_glb(3, :)' * simplified_h1234 + ...
    cord_glb(6, :)' * simplified_h56 + ...
    cord_glb(7, :)' * h(7) + ...
    cord_glb(8, :)' * h(8);

simplify_X = simplify(X);
clear X i

J = [diff(simplify_X, r), diff(simplify_X, s), diff(simplify_X, t)].';
J_inv = inv(J);

dNdr = [diff(simplified_h1234, r), diff(simplified_h56, r), diff(h(7), r), diff(h(8), r)];
dNds = [diff(simplified_h1234, s), diff(simplified_h56, s), diff(h(7), s), diff(h(8), s)];
dNdt = [diff(simplified_h1234, t), diff(simplified_h56, t), diff(h(7), t), diff(h(8), t)];
DNDx = J \ [dNdr; dNds; dNdt];

bb = [1, 4, 6;
      4, 2, 5;
      6, 5, 3]; %应变向量位置

B = zeros(max(bb, [], 'all'), 3 * size(DNDx, 2));

for n = 1:3 % u, v, w

    for i = 1:size(DNDx, 1) % x, y, z

        for j = 1:size(DNDx, 2) % 各个节点
            B(bb(n, i), 3 * (j -1) + n) = DNDx(i, j) + B(bb(n, i), 3 * (j -1) + n);
        end

    end

end

% 后续可以进行高斯点的Bij计算，subs函数即可
% 替换实际数值
% substitutedExpr = subs(expr, {x, y}, {2, 3});

function g = g(r,s,t,cord_i)
    g = G(r, cord_i(1)) * G(s, cord_i(2)) * G(t, cord_i(3));
end

function G = G(beta, beta_i)
    G = 1/2 * (1 + beta_i * beta);
end
