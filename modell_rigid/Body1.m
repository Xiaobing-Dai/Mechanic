function [M_1,F_1,Body1,Basic_Point] = Body1(q,dq,u,L)
T_q1_q = zeros(3,6);
T_q1_q(1,5) = 1;
T_dq1_dt_q = zeros(3,6);

phi = q(1);

qe = T_q1_q' * q;
dqe = T_q1_q' * dq;

m_tot = 100;   %系统质量参数
rho = 1;
k_theta_0 = eye(3);
k_theta_0(1,1) = 100;
k_theta_0(2,2) = 33.3;
k_theta_0(3,3) = 133.33;
k_r_0c = [0.666666666666667,0.866,0]';
% u = zeros(6,1);
[M_e1,F_e1,Body1] = Mass_Force_element(qe,dqe,m_tot,rho,k_theta_0,k_r_0c);


M_1 = T_q1_q * M_e1 * T_q1_q';
F_1 = -T_q1_q * (M_e1 * T_dq1_dt_q' * dq - F_e1);
% f?r jeder Body ein minus geaddet. Achtung
k_r_00 = [0,0,0]';
Basic_Point.k_r = k_r_00;
Basic_Point.r = Body1.R * k_r_00;
Basic_Point.J_dr_dqe = zeros(6,3);
Basic_Point.J_ddr_dqedt = zeros(6,3);

Body1.dqe_dq = T_q1_q;
Body1.ddqe_dqdt = T_dq1_dt_q;

% k_r_A0 = [-0.866,1.5,0]';
% Body1.k_r_A0 = k_r_A0;
% J_dr_A0_dqe = -Body1.R * skew(k_r_A0) * Body1.R' * Body1.T_dphi_2_omega * [zeros(3,3),eye(3)];
% J_dr_A0_dqe = J_dr_A0_dqe'
% J_dr_A0_dqe = [zeros(3,3),eye(3)]' * Body1.T_dphi_2_omega' * Body1.R * skew(k_r_A0) * Body1.R';
% J_dr_A0_dqe = J_dr_A0_dqe
% Point_A0 = Point_on_Body(Body1,Basic_Point,k_r_A0);
% J_dr_A0_dq = J_dr_A0_dqe' * T_q1_q';
% Body1.J_dr_A0_dq = J_dr_A0_dq';
% J_ddr_A0_dt_dq = [sin(phi),0,-cos(phi);zeros(1,3);cos(phi),0,sin(phi)] * skew(k_r_A0) * Body1.R' * Body1.T_dphi_2_omega * [zeros(3,3),eye(3)] * T_q1_q';
% Body1.J_ddr_A0_dt_dq = J_ddr_A0_dt_dq;
% 
% k_r_B0 = [0.866,0.5,0]';
% Body1.k_r_B0 = k_r_B0;
% J_dr_B0_dq = -Body1.R * skew(k_r_B0) * [[0;1;0],zeros(3,2)];
% Body1.J_dr_B0_dq = J_dr_B0_dq';
% J_ddr_B0_dt_dq = [sin(phi),0,-cos(phi);zeros(1,3);cos(phi),0,sin(phi)] * skew(k_r_B0) * [[0;1;0],zeros(3,2)];
% Body1.J_ddr_B0_dt_dq = J_ddr_B0_dt_dq;

end
