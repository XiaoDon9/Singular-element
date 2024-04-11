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

% h236 = h(2)+ h(3) +h(6);
% simplified_h236 = simplify(h236);
% simplified_h = simplify(h);

cord_loc = [[1, -1, -1, 1, 0, -1, 0, 1].', [1, 1, -1, -1, 1, 0, -1, 0].'];

for i = 1:8
    h(i) = Ni(r, s, cord_loc(i, :));
end

simplified_h = simplify(h);

% simplify(h - h2)
% % ans
% % [0, 0, 0, 0, 0, 0, 0, 0]

cord_glb = cord_loc+1;
cord_glb(5, 1) = cord_glb(5, 1) + 0.5;
cord_glb(8, 2) = cord_glb(8, 2) + 0.5;

X = cord_glb.' * simplified_h.';
simplify_X = simplify(X);
clear X i h

J = [diff(simplify_X, r), diff(simplify_X, s)];
J_inv = inv(J);

substitutedExpr = simplify(subs(J_inv, {r}, {1}));

dNdr = diff(simplified_h, r);
dNds = diff(simplified_h, s);

DNDx = J \ [dNdr; dNds];
clear dNdr dNds

bb = [1, 3;
      3, 2]; %应变向量中偏导的位置
nGs = 9; % 3*3阶的高斯积分
mGs_rs = [
          0, 0;
          0.774596669241483, -0.774596669241483;
          0.774596669241483, 0.774596669241483;
          -0.774596669241483, 0.774596669241483;
          -0.774596669241483, -0.774596669241483;
          0, 0.774596669241483;
          0, -0.774596669241483;
          0.774596669241483, 0;
          -0.774596669241483, 0;
          ];
mGs_alpha = [
             0.888888888888889, 0.888888888888889;
             0.555555555555556, 0.555555555555556;
             0.555555555555556, 0.555555555555556;
             0.555555555555556, 0.555555555555556;
             0.555555555555556, 0.555555555555556;
             0.888888888888889, 0.555555555555556;
             0.888888888888889, 0.555555555555556;
             0.555555555555556, 0.888888888888889;
             0.555555555555556, 0.888888888888889;
             ];

BIJ = zeros(max(bb, [], 'all') * nGs, 2 * size(DNDx, 2)); 
fprintf('\nBIJ是由高斯点处的应变-位移矩阵摞起来的存储矩阵\n')
for Gsnum = 1:nGs % 各个高斯点的值

    for n = 1:2 % u, v

        for i = 1:size(DNDx, 1) % x, y
            fprintf('\n第%d个高斯点处',Gsnum)
            if n == 1
                fprintf('u')
            else
                fprintf('v')
            end

            if i == 1
                fprintf('对x的导数')
            else
                fprintf('对y的导数')
            end

            for j = 1:size(DNDx, 2) % 对应节点

                if j == 1
                    fprintf('(第一个节点贡献)\n\n放入BIJ的位置为：(%d,%d)\n', ...
                        bb(n, i) + (Gsnum - 1) * 3, 2 * (j -1) + n)
                end

                DNDx_ij = double(subs(DNDx(i, j), {r, s}, {mGs_rs(Gsnum, 1), mGs_rs(Gsnum, 2)}));
                BIJ(bb(n, i) + (Gsnum - 1) * 3, 2 * (j -1) + n) = ...
                    DNDx_ij + BIJ(bb(n, i) + (Gsnum - 1) * 3, 2 * (j -1) + n);
            end

        end

    end

end

clear i j n Gsnum

E = 1;
v = 0.3;
C = E / (1 - v ^ 2) * [1, v, 0;
                       v, 1, 0;
                       0, 0, (1 - v) / 2];

K = zeros(size(BIJ, 2), size(BIJ, 2));

for i = 1:nGs
    alpha_ij = mGs_alpha(i, 1) * mGs_alpha(i, 2);
    J_ij = double(subs(J, {r, s}, {mGs_rs(i, 1), mGs_rs(i, 2)}));
    F_ij = BIJ((i - 1) * 3 + 1:i * 3, :).' * C * BIJ((i - 1) * 3 + 1:i * 3, :) * det(J_ij);
    K = K + alpha_ij * F_ij;

end

clear alpha_ij i J_ij F_ij

function Ni = Ni(r, s, cord)
    ri = cord(1);
    si = cord(2);
    Ni = ((1 + r * ri) * (1 + s * si) - (1 - r ^ 2) * (1 + s * si) - ...
        (1 - s ^ 2) * (1 + r * ri)) * ri ^ 2 * si ^ 2/4 + ...
        (1 - r ^ 2) * (1 + s * si) * (1 - ri ^ 2) * si ^ 2/2 + ...
        (1 - s ^ 2) * (1 + r * ri) * (1 - si ^ 2) * ri ^ 2/2;
end
