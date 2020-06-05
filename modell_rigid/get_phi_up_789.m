function [phi_7_3,phi_8_3,phi_9_3] = get_phi_up_789(q,dq,L,others)
phi = q(1);
s20 = -0.6;
s2 = s20 + q(3);  %%
dphi_dt = dq(1);
dphi_dq = [1,0,0]';

Body2 = others{1};
Body6 = others{2};
Point_D = others{3};
Point_H = others{4};
phi_2 = others{5};

L_8 = L.L_8;
L_9 = L.L_9;
k5_r_DF = L.k5_r_F;
k5_r_DL = L.k5_r_L;
k7_r_HF = L.k7_r_H;
k8_r_HL = [L_8 + L_9 + s2;0;0];

r_H = Point_H.r;
r_D = Point_D.r;
%%
[phi_3,f] = fmincon(@(phi)get_phi_3_func(phi,s2,L,others),[0,0],[],[],[],[],[-pi,0],[2 * pi,pi]);
phi_27_3 = phi_3(1);
phi_28_3 = phi_3(2);
phi_7_3.phi_27_3 = phi_27_3;
phi_7_3.phi_7_3 = phi_2.phi_2_3 + phi_27_3;
phi_8_3.phi_28_3 = phi_28_3;
phi_8_3.phi_8_3 = phi_2.phi_2_3 + phi_28_3;
%%
clear A B C D dD_dq BB E F G H J;
R_2 = Body2.R;
R_28 = [cos(phi_28_3),-sin(phi_28_3),0;sin(phi_28_3),cos(phi_28_3),0;0,0,1];
T_28_phi3 = [-sin(phi_28_3),-cos(phi_28_3),0;cos(phi_28_3),-sin(phi_28_3),0;0,0,0];

dphi_2_3_dq = [0;phi_2.dphi_2_3_ds1;0];
phi_2_3 = phi_2.phi_2_3;

R_2_2 = [cos(phi),0,sin(phi);0,0,0;-sin(phi),0,cos(phi)];
R_3_2 = [cos(phi_2_3),-sin(phi_2_3),0;sin(phi_2_3),cos(phi_2_3),0;0,0,1];
T_2_phi2 = [-sin(phi),0,cos(phi);0,0,0;-cos(phi),0,-sin(phi)];
T_2_phi3 = [-sin(phi_2_3),-cos(phi_2_3),0;cos(phi_2_3),-sin(phi_2_3),0;0,0,0];

dqe2_dq = Body2.dqe_dq;
dqe6_dq = Body6.dqe_dq;
dr_D_dqe2 = Point_D.J_dr_dqe;
dr_H_dqe6 = Point_H.J_dr_dqe;

dr_H_dq = dqe6_dq * dr_H_dqe6;
dr_D_dq = dqe2_dq * dr_D_dqe2;

A = r_H - r_D + R_2 * R_28 * k8_r_HL;
B = T_2_phi2 * R_3_2 * R_28 * k8_r_HL * dphi_dq' + R_2_2 * T_2_phi3 * R_28 * k8_r_HL * dphi_2_3_dq';
C = -A' * (R_2 * T_28_phi3 * k8_r_HL);

dphi_28_3_dq = 1/C * A' * (dr_H_dq' - dr_D_dq' + B);
dphi_28_3_dq = dphi_28_3_dq';
dphi_8_3_dq = phi_2.dphi_2_3_dq + dphi_28_3_dq;
phi_8_3.dphi_8_3_dq = dphi_8_3_dq;
%%
dphi_2_3_dt = dphi_2_3_dq' * dq;
dphi_8_3_dt = dphi_8_3_dq' * dq;

dR_2_dt = R_2_2 * T_2_phi3 * dphi_2_3_dt + T_2_phi2 * R_3_2 * dphi_dt;
dR_28_dt = T_28_phi3 * dphi_8_3_dt;

dT_2_phi2_dt = [-cos(phi),0,-sin(phi);0,0,0;sin(phi),0,-cos(phi)] * dphi_dt;
dT_2_phi3_dt = [-cos(phi_2_3),sin(phi_2_3),0;-sin(phi_2_3),-cos(phi_2_3),0;0,0,0] * dphi_2_3_dt;

dr_H_dt = dr_H_dq' * dq;
dr_D_dt = dr_D_dq' * dq;
ddr_H_dq = Body6.dqe_dq * Point_H.J_ddr_dqedt + Body6.ddqe_dqdt * Point_H.J_dr_dqe;
ddr_D_dq = Body2.dqe_dq * Point_D.J_ddr_dqedt + Body2.ddqe_dqdt * Point_D.J_dr_dqe;

D = dr_H_dt - dr_D_dt + dR_2_dt * R_28 * k8_r_HL + R_2 * dR_28_dt * k8_r_HL;
dD_dq = dr_H_dq' - dr_D_dq' + (R_2_2 * T_2_phi3 * R_28 * k8_r_HL * dphi_2_3_dq' + T_2_phi2 * R_3_2 * R_28 * k8_r_HL * dphi_dq') + ...
    R_2 * T_28_phi3 * k8_r_HL * dphi_8_3_dq';
dD_dq = dD_dq';

BB = T_2_phi2 * R_3_2 * dR_28_dt * k8_r_HL * dphi_dq' + R_2_2 * T_2_phi3 * dR_2_dt * k8_r_HL * dphi_2_3_dq';
E = dT_2_phi2_dt * R_3_2 * R_28 * k8_r_HL * dphi_dq' + 2 * T_2_phi2 * T_2_phi3 * R_28 * k8_r_HL * dphi_dq' * dphi_2_3_dt + ...
    R_2_2 * dT_2_phi3_dt * R_28 * k8_r_HL * dphi_2_3_dq' + R_2_2 * T_2_phi3 * R_28 * k8_r_HL * dphi_2_3_dq';
F = T_2_phi2 * R_3_2 * R_28 * k8_r_HL * dphi_dq' + R_2_2 * T_2_phi3 * R_28 * k8_r_HL * dphi_2_3_dq';
G = -A' * R_2 * T_28_phi3 * k8_r_HL;
H = A' * (ddr_H_dq' - ddr_D_dq' + E + 2 * BB) + D' * dD_dq';
H = H';
J = A' * (dr_H_dq' - dr_D_dq' + F);
J = J';

ddphi_28_3_dqdt = 1/G * H';
ddphi_28_3_dqdt = ddphi_28_3_dqdt';
ddphi_8_3_dqdt = phi_2.ddphi_2_3_dqdt + ddphi_28_3_dqdt;
phi_8_3.ddphi_8_3_dqdt = ddphi_8_3_dqdt;
clear A B C D dD_dq BB E F G H J;
%%
clear A dA_dq dA_dt B dB_dq C E;
R_2_2 = [cos(phi),0,sin(phi);0,1,0;-sin(phi),0,cos(phi)];
R_2_3 = [cos(phi_2_3),-sin(phi_2_3),0;sin(phi_2_3),cos(phi_2_3),0;0,0,1];

T_2_phi2 = [-sin(phi),0,cos(phi);0,0,0;-cos(phi),0,-sin(phi)];
T_2_phi3 = [-sin(phi_2_3),-cos(phi_2_3),0;cos(phi_2_3),-sin(phi_2_3),0;0,0,0];

dR_2_dt = Body2.dR_dt;

R_27 = [cos(phi_27_3),-sin(phi_27_3),0;sin(phi_27_3),cos(phi_27_3),0;0,0,1];
T_27_phi3 = [-sin(phi_27_3),-cos(phi_27_3),0;cos(phi_27_3),-sin(phi_27_3),0;0,0,0];



A = r_H - r_D + R_2 * R_27 * k7_r_HF;
dA_dq = (-Point_D.J_dr_dqe' * Body2.dqe_dq' + Point_H.J_dr_dqe' * Body6.dqe_dq' + ...
    T_2_phi2 * R_2_3 * R_27 * k7_r_HF * dphi_dq' + ...
    R_2_2 * T_2_phi3 * R_27 * k7_r_HF * phi_2.dphi_2_3_dq')';
B = R_2 * T_27_phi3 * k7_r_HF;


%
dphi_27_3_dq = - 1 / (A' * B) * dA_dq * A;
dphi_7_3_dq = phi_2.dphi_2_3_dq + dphi_27_3_dq;
phi_7_3.dphi_7_3_dq = dphi_7_3_dq;
%
dT_27_phi3_dphi_27_3 = [-cos(phi_27_3),sin(phi_27_3),0;-sin(phi_27_3),-cos(phi_27_3),0;0,0,0];
dT_27_phi3_dt = dT_27_phi3_dphi_27_3 * (dphi_27_3_dq' * dq);

dA_dt = dA_dq' * dq;
dB_dq = T_2_phi2 * R_2_3 * T_27_phi3 * k7_r_HF * dphi_dq' ...
      + R_2_2 * T_2_phi3 * T_27_phi3 * k7_r_HF * phi_2.dphi_2_3_dq' ...
      + R_2 * dT_27_phi3_dphi_27_3 * k7_r_HF * dphi_27_3_dq';
dB_dq = dB_dq';
% dB_dt = dB_dq' * dq;

dR_2_2_dt = T_2_phi2 * (dphi_dq' * dq);
dR_3_2_dt = T_2_phi3 * (phi_2.dphi_2_3_dq' * dq);

C = (dT_2_phi2_dt * R_3_2 + T_2_phi2 * dR_3_2_dt) * R_27 * k7_r_HF * dphi_dq'+ ...
    (dR_2_2_dt * T_2_phi3 + R_2_2 * dT_2_phi3_dt) * R_27 * k7_r_HF * phi_2.dphi_2_3_dq' + ...
    R_2_2 * T_2_phi3 * R_27 * k7_r_HF * phi_2.ddphi_2_3_dqdt';
E = ddr_H_dq' - ddr_D_dq' + C + dR_2_dt * T_27_phi3 * k7_r_HF * dphi_27_3_dq';


ddphi_27_3_dqdt = -1 / (A' * B) * (dA_dt' * dA_dq' + A' * E + B' * dA_dq' + A' * dB_dq');
ddphi_27_3_dqdt = ddphi_27_3_dqdt';
ddphi_7_3_dqdt = phi_2.ddphi_2_3_dqdt + ddphi_27_3_dqdt;
phi_7_3.ddphi_7_3_dqdt = ddphi_7_3_dqdt;

%%
phi_9_3.phi_29_3 = phi_28_3;
phi_9_3.phi_9_3 = phi_2.phi_2_3 + phi_28_3;
phi_9_3.dphi_9_3_dq = dphi_8_3_dq;
phi_9_3.ddphi_9_3_dqdt = ddphi_8_3_dqdt;

end

function tolerance = get_phi_3_func(phi_3,s2,L,others)
phi_27_3 = phi_3(1);
phi_28_3 = phi_3(2);

Body2 = others{1};
Point_D = others{3};
Point_H = others{4};

r_H = Point_H.r;
r_D = Point_D.r;
R_2 = Body2.R;

L_8 = L.L_8;
L_9 = L.L_9;
k5_r_DF = L.k5_r_F;
k5_r_DL = L.k5_r_L;
k7_r_HF = L.k7_r_H;
k8_r_HL = [L_8 + L_9 + s2;0;0];

R_27 = [cos(phi_27_3),-sin(phi_27_3),0;sin(phi_27_3),cos(phi_27_3),0;0,0,1];
R_28 = [cos(phi_28_3),-sin(phi_28_3),0;sin(phi_28_3),cos(phi_28_3),0;0,0,1];

RHS1 = norm(k5_r_DF);
RHS2 = norm(k5_r_DL);

LHS1 = norm(r_H - r_D + R_2 * R_27 * k7_r_HF);
LHS2 = norm(r_H - r_D + R_2 * R_28 * k8_r_HL);

tolerance = (RHS1 - LHS1) ^ 2 + (RHS2 - LHS2) ^ 2;

end