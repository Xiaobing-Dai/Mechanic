function [M_3,F_3,Body3] = Body3(q,dq,u,others,L)
Body2 = others{1};
Body4 = others{2};
Point_A = others{3};
J_dr_A_dq = Body2.dqe_dq * Point_A.J_dr_dqe;
J_ddr_A_dq = Body2.dqe_dq * Point_A.J_ddr_dqedt + Body2.ddqe_dqdt * Point_A.J_dr_dqe;
dphi_4_3_ds1 = Body4.dphi_4_3_ds1;
ddphi_4_3_dt_ds1 = Body4.ddphi_4_3_dt_ds1;

qe = zeros(6,1);
qe(1:3,1) = Point_A.r;
qe(4:6,1) = Body4.phi;

T_q3_q = [J_dr_A_dq';zeros(1,3);1,0,0;0,dphi_4_3_ds1,0]';
dqe = T_q3_q' * dq;

T_dq3_dt_q = [J_ddr_A_dq';zeros(1,3);1,0,0;0,ddphi_4_3_dt_ds1,0]';

m_tot = 10;
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 0.0061;
k_theta_0(2,2) = 0.8364;
k_theta_0(3,3) = 0.8364;
k_r_0c = [-0.5,0,0]';
% u = zeros(6,1);
[M_e3,F_e3,Body3] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);
%% Õ‚¡¶
phi_1 = qe(4);
phi_2 = qe(5);
phi_3 = qe(6);
R_1 = [1,0,0;0,cos(phi_1),-sin(phi_1);0,sin(phi_1),cos(phi_1)];
R_2 = [cos(phi_2),0,sin(phi_2);0,1,0;-sin(phi_2),0,cos(phi_2)];
R_3 = [cos(phi_3),-sin(phi_3),0;sin(phi_3),cos(phi_3),0;0,0,1];
R = R_2 * R_3 * R_1;
F = u(2) * R * [1;0;0];
F_k = [eye(3)*F;skew(Point_A.r + R * L.k3_AB)*F];
F_e3 = F_e3 + F_k;
%%
M_3 = T_q3_q * M_e3 * T_q3_q';
F_3 = -T_q3_q * (M_e3 * T_dq3_dt_q' * dq - F_e3);

end



