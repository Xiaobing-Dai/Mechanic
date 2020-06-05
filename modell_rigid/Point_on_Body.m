function Point = Point_on_Body(Body,Basic_point,k_r_0A)
T_phiq = [zeros(3),eye(3)];
J_dr_0_dqe = Basic_point.J_dr_dqe;
J_ddr_0_dqedt = Basic_point.J_ddr_dqedt;
T_dphi_2_omega = Body.T_dphi_2_omega;
R = Body.R;
omega = Body.omega;

J_dr_A_dqe = J_dr_0_dqe + T_phiq' * T_dphi_2_omega' * R * skew(k_r_0A) * R';
J_ddr_A_dqedt = J_ddr_0_dqedt - T_phiq' * T_dphi_2_omega' * R * skew(k_r_0A) * skew(omega) * R';
%- T_phiq'  ein minus geaddet
Point.J_dr_dqe = J_dr_A_dqe;
Point.J_ddr_dqedt = J_ddr_A_dqedt;
Point.k_r = k_r_0A;
Point.r = Basic_point.r + R * k_r_0A;
end