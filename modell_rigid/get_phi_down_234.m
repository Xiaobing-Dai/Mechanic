function phi_3 = get_phi_down_234(q,L)
%%
rA0 = L.k1_r_0A;
x0 = rA0(1);y0 = rA0(2);z0 = rA0(3);

rB0 = L.k1_r_0B;
x1 = rB0(1);y1 = rB0(2);z1 = rB0(3);

L1 = norm(rA0 - rB0);
L2 = abs(L.k2_r_A(1));
L3 = L.L_3;
L4 = L.L_4;
s1 = q(2);

Q = ((L3+L4+s1)^2 - L1^2 - L2^2) / (2*L2);

A = Q + (x0 - x1);
B = -2 * (y0 - y1);
C = Q - (x0 - x1);

p1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A);
p2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A);

phi2 = atan(p2) * 2;


%%
rA = rA0 + [L2*cos(phi2);L2*sin(phi2);0];
delta = rA - rB0;
deltax = delta(1);
deltay = delta(2);
if rA(1) >= rB0(1)
	phi4 = atan(deltay / deltax);
else
	phi4 = pi + atan(deltay / deltax);
end
%%
phi_3 = [phi2,phi4];
end