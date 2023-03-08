format long;
r = 4;
[p, pot, x, y] = Gibanje_v_krogu(r, 4, 0.1, 100, pi);
verjetnost = p
ezplot(@(x,y) (x).^2 + (y).^2 -r^2);
hold on;
plot(pot(1,:), pot(2,:), 'MarkerFaceColor', 'k','LineWidth', 1.5);
hold on;
%plot(x,y, 'r*');
