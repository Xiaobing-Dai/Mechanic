function [M_6,F_6,Body6] = Body6(q,dq,u,others,L)
Body2 = others{1};
Point_E = others{2};
phi_6_3 = others{3};

dphi_2_3_ds1 = Body2.dphi_2_3_ds1;
ddphi_2_3_ds1dt = Body2.ddphi_2_3_dt_ds1;
dphi_26_3_ds2 = phi_6_3.dphi_6_3_ds2;
ddphi_26_3_ds2dt = phi_6_3.ddphi_6_3_ds2dt;

dr_E_dqe = Body2.dqe_dq * Point_E.J_dr_dqe;
dqe_6_dq = [dr_E_dqe';zeros(1,3);1,0,0;0,dphi_2_3_ds1,dphi_26_3_ds2]';

ddr_E_dqdt = Body2.dqe_dq * Point_E.J_ddr_dqedt + Body2.ddqe_dqdt * Point_E.J_dr_dqe;
ddqe_6_dqdt = [ddr_E_dqdt';zeros(1,3);1,0,0;0,ddphi_2_3_ds1dt,ddphi_26_3_ds2dt]';

phi_6 = Body2.phi + [0;0;phi_6_3.phi_26_3];
qe = zeros(6,1);
qe(1:3,1) = Point_E.r;
qe(4:6,1) = phi_6;
dqe = dqe_6_dq' * dq;

m_tot = 10;
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 0.0125;
k_theta_0(2,2) = 0.43;
k_theta_0(3,3) = 0.43;
k_r_0c = [0.875,0,0]';
% u = zeros(6,1);
[M_e6,F_e6,Body6] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);

M_6 = dqe_6_dq * M_e6 * dqe_6_dq';
F_6 = -dqe_6_dq * (M_e6 * ddqe_6_dqdt' * dq - F_e6);

Body6.dqe_dq = dqe_6_dq;
Body6.ddqe_dqdt = ddqe_6_dqdt;
end