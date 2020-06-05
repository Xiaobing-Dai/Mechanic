function [M_7,F_7,Body7] = Body7(q,dq,u,others,L)
phi = q(1);

Body6 = others{1};
Point_H = others{2};
phi_7_3 = others{3};

phi_7 = [0;phi;phi_7_3.phi_7_3];

qe = zeros(6,1);
qe(1:3,1) = Point_H.r;
qe(4:6,1) = phi_7;

dr_H_dq = Body6.dqe_dq * Point_H.J_dr_dqe;
dqe_7_dq = [dr_H_dq';zeros(1,3);1,0,0;phi_7_3.dphi_7_3_dq']';
dqe = dqe_7_dq' * dq;

ddr_H_dqdt = Body6.dqe_dq * Point_H.J_ddr_dqedt + Body6.ddqe_dqdt * Point_H.J_dr_dqe;
ddqe_7_dqdt = [ddr_H_dqdt';zeros(1,3);1,0,0;phi_7_3.ddphi_7_3_dqdt']';

m_tot = 5;
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 0.0125;
k_theta_0(2,2) = 0.43;
k_theta_0(3,3) = 0.43;
k_r_0c = [0.875,0,0]';
% u = zeros(6,1);
[M_e7,F_e7,Body7] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);

M_7 = dqe_7_dq * M_e7 * dqe_7_dq';
F_7 = -dqe_7_dq * (M_e7 * ddqe_7_dqdt' * dq - F_e7);
end
