function [M_9,F_9,Body9] = Body9(q,dq,u,others,L)
phi = q(1);

Body5 = others{1};
Point_L = others{2};
phi_9_3 = others{3};

phi_9 = [0;phi;phi_9_3.phi_9_3];

qe = zeros(6,1);
qe(1:3,1) = Point_L.r;
qe(4:6,1) = phi_9;

dr_H_dq = Body5.dqe_dq * Point_L.J_dr_dqe;
dqe_9_dq = [dr_H_dq';zeros(1,3);1,0,0;phi_9_3.dphi_9_3_dq']';
dqe = dqe_9_dq' * dq;

ddr_H_dqdt = Body5.dqe_dq * Point_L.J_ddr_dqedt + Body5.ddqe_dqdt * Point_L.J_dr_dqe;
ddqe_9_dqdt = [ddr_H_dqdt';zeros(1,3);1,0,0;phi_9_3.ddphi_9_3_dqdt']';

m_tot = 10;
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 0.0061;
k_theta_0(2,2) = 0.8364;
k_theta_0(3,3) = 0.8364;
k_r_0c = [-0.5,0,0]';
% u = zeros(6,1);
[M_e9,F_e9,Body9] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);
%%
phi_1 = qe(4);
phi_2 = qe(5);
phi_3 = qe(6);
R_1 = [1,0,0;0,cos(phi_1),-sin(phi_1);0,sin(phi_1),cos(phi_1)];
R_2 = [cos(phi_2),0,sin(phi_2);0,1,0;-sin(phi_2),0,cos(phi_2)];
R_3 = [cos(phi_3),-sin(phi_3),0;sin(phi_3),cos(phi_3),0;0,0,1];
R = R_2 * R_3 * R_1;
F = u(3) * R * [1;0;0];
F_k = [eye(3)*F;skew(Point_L.r + R * L.k9_LG)*F];
F_e9 = F_e9 + F_k;
%%
M_9 = dqe_9_dq * M_e9 * dqe_9_dq';
F_9 = -dqe_9_dq * (M_e9 * ddqe_9_dqdt' * dq - F_e9);
end
