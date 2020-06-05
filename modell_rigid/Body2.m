function [M_2,F_2,Body2,phi_2_3] = Body2(q,dq,u,others,L)
phi = q(1);
s1 = q(2);
ds1 = dq(2);
Body1 = others{1};
Point_A0 = others{2};
Point_B0 = others{3};

L_3 = L.L_3;
L_4 = L.L_4;
L_s1 = L_3 + L_4 + s1;
L_2a = L.L_A0A;

k_r_B0A0 = Point_A0.k_r - Point_B0.k_r;
L_1 = norm(k_r_B0A0);
C = L_s1 ^ 2 - L_1 ^ 2 - L_2a ^ 2;
A = 2 * k_r_B0A0(1) * L_2a;
B = 2 * k_r_B0A0(2) * L_2a;
phi_2 = [0,phi,0]';
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
phi_3 = get_phi_down_234(q,L);
% [phi_3,f] = fminunc(@(phi)sum([k_r_B0A0(1) + L_2a * cos(phi(1)) - L_s1 * cos(phi(2));...
%     k_r_B0A0(2) + L_2a * sin(phi(1)) - L_s1 * sin(phi(2))] .^ 2),[0,0],options);

phi_3_12 = phi_3(1);
% phi_3_14 = phi_3(2);
phi_2(3) = phi_3_12;
qe = zeros(6,1);
qe(1:3,1) = Point_A0.r;
qe(4:6,1) = phi_2;

J_dr_A0_dq = Body1.dqe_dq * Point_A0.J_dr_dqe;
J_ddr_A0_dq = Body1.dqe_dq * Point_A0.J_ddr_dqedt + Body1.ddqe_dqdt * Point_A0.J_dr_dqe;
dphi_2_3_ds1 = -L_s1 / (L_2a * ((k_r_B0A0(1) * sin(phi_2(3))) - (k_r_B0A0(2) * cos(phi_2(3)))));
ddphi_2_3_dt_ds1 = -ds1 / L_2a / (k_r_B0A0(1) * sin(phi_2(3)) - k_r_B0A0(2) * cos(phi_2(3))) - ...
    ds1 * (L_s1 ^ 2) * (k_r_B0A0(1) * cos(phi_2(3)) + k_r_B0A0(2) * sin(phi_2(3))) / (L_2a ^ 2) / ((k_r_B0A0(1) * sin(phi_2(3)) - k_r_B0A0(2) * cos(phi_2(3))) ^ 3);
T_q2_q = [J_dr_A0_dq';zeros(1,3);1,0,0;0,dphi_2_3_ds1,0]';
dqe = T_q2_q' * dq;

T_dq2_dt_q = [J_ddr_A0_dq;zeros(1,3);1,0,0;0,ddphi_2_3_dt_ds1,0]';

m_tot = 50;   %
rho = 1;
k_theta_0 = eye(3) ;  %转动惯量
k_theta_0(1,1) = 0.225;    
k_theta_0(2,2) = 5.1125;
k_theta_0(3,3) = 5.1125;
k_r_0c = [0.866,0,0]';
% k_r_0c = [13'/2,0,0]';    %质心位置
% u = zeros(6,1);
[M_e2,F_e2,Body2] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);

M_2 = T_q2_q * M_e2 * T_q2_q';
F_2 = -T_q2_q * (M_e2 * T_dq2_dt_q' * dq - F_e2);

Body2.dqe_dq = T_q2_q;
Body2.ddqe_dqdt = T_dq2_dt_q;
Body2.dphi_2_3_ds1 = dphi_2_3_ds1;
Body2.ddphi_2_3_dt_ds1 = ddphi_2_3_dt_ds1;
Body2.phi = phi_2;

phi_2_3.phi_2_3 = phi_3_12;
phi_2_3.dphi_2_3_ds1 = dphi_2_3_ds1;
phi_2_3.dphi_2_3_dq = [0;dphi_2_3_ds1;0];
phi_2_3.ddphi_2_3_dqdt = [0;ddphi_2_3_dt_ds1;0];
end










