function [M_5,F_5,Body5] = Body5(q,dq,u,others,L)
Body2 = others{1};
Point_D = others{2};
phi_5_3 = others{3};

dphi_2_3_ds1 = Body2.dphi_2_3_ds1;
ddphi_2_3_ds1dt = Body2.ddphi_2_3_dt_ds1;
dphi_25_3_ds2 = phi_5_3.dphi_5_3_ds2;
ddphi_25_3_ds2dt = phi_5_3.ddphi_5_3_ds2dt;

dr_D_dq = Body2.dqe_dq * Point_D.J_dr_dqe;
dqe_5_dq = [dr_D_dq';zeros(1,3);1,0,0;0,dphi_2_3_ds1,dphi_25_3_ds2]';

ddr_D_dqdt = Body2.dqe_dq * Point_D.J_ddr_dqedt + Body2.ddqe_dqdt * Point_D.J_dr_dqe;
ddqe_5_dqdt = [ddr_D_dqdt';zeros(1,3);1,0,0;0,ddphi_2_3_ds1dt,ddphi_25_3_ds2dt]';

phi_5 = Body2.phi + [0;0;phi_5_3.phi_25_3];
qe = zeros(6,1);
qe(1:3,1) = Point_D.r;
qe(4:6,1) = phi_5;
dqe = dqe_5_dq' * dq;

m_tot = 50;
rho = 1;
k_theta_0 = eye(3) ;
k_theta_0(1,1) = 11.25;
k_theta_0(2,2) = 8339;
k_theta_0(3,3) = 8339;
k_r_0c = [0.5,0,0]';
% u = zeros(6,1);
[M_e5,F_e5,Body5] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);

M_5 = dqe_5_dq * M_e5 * dqe_5_dq';
F_5 = -dqe_5_dq * (M_e5 * ddqe_5_dqdt' * dq - F_e5);

Body5.dqe_dq = dqe_5_dq;
Body5.ddqe_dqdt = ddqe_5_dqdt;
end