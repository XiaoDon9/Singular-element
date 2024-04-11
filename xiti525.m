clear
clc

syms r s t

h(8) = 1/2 * (1 - s^2) * (1 + r);
h(7) = 1/2 * (1 - r^2) * (1 - s);
h(6) = 1/2 * (1 - s^2) * (1 - r);
h(5) = 1/2 * (1 - r^2) * (1 + s);

h(4) = 1/4 * (1 + r) * (1 - s) - 1/2 * h(7) - 1/2 * h(8);
h(3) = 1/4 * (1 - r) * (1 - s) - 1/2 * h(6) - 1/2 * h(7); 
h(2) = 1/4 * (1 - r) * (1 + s) - 1/2 * h(5) - 1/2 * h(6); 
h(1) = 1/4 * (1 + r) * (1 + s) - 1/2 * h(5) - 1/2 * h(8); 

h125 = h(1) + h(2) + h(5);
simplified_h125 = simplify(h125);

delta_h = 1/8 * (1 - r^2) * (1 - s^2);
h(3) = h(3) + delta_h;
h(4) = h(4) + delta_h;
h(7) = h(7) + delta_h;
simplified_h = simplify(h);

disp(simplified_h125)
disp(simplified_h([3,4,6,7,8]).')

latex(simplified_h125)
latex(simplified_h([3,4,6,7,8]).')

symdisp(simplified_h125)
symdisp(simplified_h([3,4,6,7,8]).')