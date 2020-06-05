function [M_4,F_4,Body4] = Body4(q,dq,u,others,L)
Body1 = others{1};
Point_A0 = others{2};
Point_B0 = others{3};
phi = q(1);
s1 = q(2);
ds1 = dq(2);

L_3 = L.L_3;
L_4 = L.L_4;
L_s1 = L_3 + L_4 + s1;
L_2a = L.L_A0A;

k_r_B0A0 = Point_A0.k_r - Point_B0.k_r;
L_1 = norm(k_r_B0A0);

phi_4 = [0,phi,0]';
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
% [phi_3,f] = fminunc(@(phi)sum([k_r_B0A0(1) + L_2a * cos(phi(1)) - L_s1 * cos(phi(2));...
%     k_r_B0A0(2) + L_2a * sin(phi(1)) - L_s1 * sin(phi(2))] .^ 2),[0,0],options);
phi_3 =  get_phi_down_234(q,L);
phi_3_14 = phi_3(2);
phi_4(3) = phi_3_14;

qe = zeros(6,1);
qe(1:3,1) = Point_B0.r;
qe(4:6,1) = phi_4;

J_dr_B0_dq = Body1.dqe_dq * Point_B0.J_dr_dqe;
J_ddr_B0_dq = Body1.dqe_dq * Point_B0.J_ddr_dqedt + Body1.ddqe_dqdt * Point_B0.J_dr_dqe;

A = L_s1 - k_r_B0A0(1) * cos(phi_4(3)) - k_r_B0A0(2) * sin(phi_4(3));
B = L_s1 * (k_r_B0A0(2) * cos(phi_4(3)) - k_r_B0A0(1) * sin(phi_4(3)));
dphi_4_3_ds1 = A / B;
T_q4_q = [J_dr_B0_dq';zeros(1,3);1,0,0;0,dphi_4_3_ds1,0]';
dqe = T_q4_q' * dq;

dA = 1 + (k_r_B0A0(1) * sin(phi_4(3)) - k_r_B0A0(2) * cos(phi_4(3))) * dphi_4_3_ds1;
dB = k_r_B0A0(2) * cos(phi_4(3)) - k_r_B0A0(1) * sin(phi_4(3)) - ...
    L_s1 * (k_r_B0A0(1) * cos(phi_4(3)) + k_r_B0A0(2) * sin(phi_4(3))) * dphi_4_3_ds1;
ddphi_4_3_dt_ds1 = ds1 / (B ^ 2) * (dA * B - A * dB);
T_dq4_dt_q = [J_ddr_B0_dq';zeros(1,3);1,0,0;0,ddphi_4_3_dt_ds1,0]';


m_tot = 10;
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 0.0125;
k_theta_0(2,2) = 0.8396;
k_theta_0(3,3) = 0.8396;
k_r_0c = [0.5,0,0]';
% u = zeros(6,1);
[M_e4,F_e4,Body4] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);
%%
phi_1 = qe(4);
phi_2 = qe(5);
phi_3 = qe(6);
R_1 = [1,0,0;0,cos(phi_1),-sin(phi_1);0,sin(phi_1),cos(phi_1)];
R_2 = [cos(phi_2),0,sin(phi_2);0,1,0;-sin(phi_2),0,cos(phi_2)];
R_3 = [cos(phi_3),-sin(phi_3),0;sin(phi_3),cos(phi_3),0;0,0,1];
R = R_2 * R_3 * R_1;
F = u(2) * R * [-1;0;0];
F_k = [eye(3)*F;skew(Point_B0.r + R * L.k3_B0BB)*F];
F_e4 = F_e4 + F_k;
%%
M_4 = T_q4_q * M_e4 * T_q4_q';
F_4 = -T_q4_q * (M_e4 * T_dq4_dt_q' * dq - F_e4);

Body4.dphi_4_3_ds1 = dphi_4_3_ds1;
Body4.ddphi_4_3_dt_ds1 = ddphi_4_3_dt_ds1;
Body4.phi = phi_4;
end