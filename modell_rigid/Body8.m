function [M_8,F_8,Body8] = Body8(q,dq,u,others,L)
phi = q(1);

Body6 = others{1};
Point_H = others{2};
phi_8_3 = others{3};

phi_8 = [0;phi;phi_8_3.phi_8_3];

qe = zeros(6,1);
qe(1:3,1) = Point_H.r;
qe(4:6,1) = phi_8;

dr_H_dq = Body6.dqe_dq * Point_H.J_dr_dqe;
dqe_8_dq = [dr_H_dq';zeros(1,3);1,0,0;phi_8_3.dphi_8_3_dq']';
dqe = dqe_8_dq' * dq;

ddr_H_dqdt = Body6.dqe_dq * Point_H.J_ddr_dqedt + Body6.ddqe_dqdt * Point_H.J_dr_dqe;
ddqe_8_dqdt = [ddr_H_dqdt';zeros(1,3);1,0,0;phi_8_3.ddphi_8_3_dqdt']';

m_tot = 10;
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 0.0125;
k_theta_0(2,2) = 0.8396;
k_theta_0(3,3) = 0.8396;
k_r_0c = [0.5,0,0]';
% u = zeros(6,1);
[M_e8,F_e8,Body8] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);
%%
phi_1 = qe(4);
phi_2 = qe(5);
phi_3 = qe(6);
R_1 = [1,0,0;0,cos(phi_1),-sin(phi_1);0,sin(phi_1),cos(phi_1)];
R_2 = [cos(phi_2),0,sin(phi_2);0,1,0;-sin(phi_2),0,cos(phi_2)];
R_3 = [cos(phi_3),-sin(phi_3),0;sin(phi_3),cos(phi_3),0;0,0,1];
R = R_2 * R_3 * R_1;
F = u(3) * R * [-1;0;0];
F_k = [eye(3)*F;skew(Point_H.r + R * L.k8_HGG)*F];
F_e8 = F_e8 + F_k;
%%
M_8 = dqe_8_dq * M_e8 * dqe_8_dq';
F_8 = -dqe_8_dq * (M_e8 * ddqe_8_dqdt' * dq - F_e8);
end
