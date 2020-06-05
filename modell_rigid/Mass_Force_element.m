function [M_e,F_e,others] = Mass_Force_element(q,dq,m_tot,rho,k_theta_0,k_r_0c)
T_rq = [eye(3),zeros(3)];
T_phiq = [zeros(3),eye(3)];
% eye ,vorher ones
phi_1 = q(4);
phi_2 = q(5);
phi_3 = q(6);
phi = T_phiq * q;
dphi = T_phiq * dq;

R_1 = [1,0,0;0,cos(phi_1),-sin(phi_1);0,sin(phi_1),cos(phi_1)];
R_2 = [cos(phi_2),0,sin(phi_2);0,1,0;-sin(phi_2),0,cos(phi_2)];
R_3 = [cos(phi_3),-sin(phi_3),0;sin(phi_3),cos(phi_3),0;0,0,1];
R = R_2 * R_3 * R_1;

dR_1_dt = [0,0,0;0,-sin(phi_1),-cos(phi_1);0,cos(phi_1),-sin(phi_1)] * dphi(1);
dR_2_dt = [-sin(phi_2),0,cos(phi_2);0,0,0;-cos(phi_2),0,-sin(phi_2)] * dphi(2);
dR_3_dt = [-sin(phi_3),-cos(phi_3),0;cos(phi_3),-sin(phi_3),0;0,0,0] * dphi(3);
dR_dt = R_2 * R_3 * dR_1_dt + dR_2_dt * R_3 * R_1 + R_2 * dR_3_dt * R_1;

theta_0 = R * k_theta_0 * R';
r_0c = R * k_r_0c;

T_dphi_2_omega = [zeros(1,3);0,1,0;zeros(1,3)] + R_2 * [zeros(2,3);0,0,1] + R_2 * R_3 * [1,0,0;zeros(2,3)];
omega = T_dphi_2_omega * dphi;


r_0c_skew = skew(r_0c);
omega_skew = skew(omega);
%%
%M_e,F_ine_e
M_e = [ m_tot * eye(3),-m_tot * r_0c_skew * T_dphi_2_omega;
        -m_tot * T_dphi_2_omega' * r_0c_skew',T_dphi_2_omega' * theta_0 * T_dphi_2_omega];
F_ine_e = (-m_tot * T_rq' * omega_skew * r_0c_skew + T_phiq' * T_dphi_2_omega' * omega_skew * theta_0) * T_dphi_2_omega * T_phiq * dq;
%%
g = 9.8;
F_ext_e_g = -(m_tot * T_rq' + m_tot / rho * T_phiq' * T_dphi_2_omega' * r_0c_skew) * [0,g,0]';
% vorher (0,0,g)
F_ext_e = F_ext_e_g;
%%
F_e = F_ext_e - F_ine_e;
% Achtung,vor F_ext_e ein minus geaddet
%%
others.R_1 = R_1;
others.R_2 = R_2;
others.R_3 = R_3;
others.R = R;
others.dR_dt = dR_dt;
others.T_dphi_2_omega = T_dphi_2_omega;
others.omega = omega;
end