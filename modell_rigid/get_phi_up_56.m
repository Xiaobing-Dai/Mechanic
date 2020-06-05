function [phi_5_3,phi_6_3] = get_phi_up_56(q,dq,L)
s20 = -0.6; %%
s2 = s20+ q(3);   %%%
ds2 = dq(3);

k2_r_D = L.k2_r_D;
k2_r_E = L.k2_r_E;
k5_r_DF = L.k5_r_F;
k5_r_DL = L.k5_r_L;
k6_r_EH = L.k6_r_H;

[phi_3,f] = fmincon(@(phi)get_phi_3_func(phi,s2,L),[pi/2,0],[],[],[],[],[0,0],[2 * pi,2 * pi]); %%

phi_25_3 = phi_3(1);
phi_26_3 = phi_3(2);

R_25 = [cos(phi_25_3),-sin(phi_25_3),0;sin(phi_25_3),cos(phi_25_3),0;0,0,1];
R_26 = [cos(phi_26_3),-sin(phi_26_3),0;sin(phi_26_3),cos(phi_26_3),0;0,0,1];

A0 = k2_r_D - k2_r_E + R_25 * k5_r_DF - R_26 * k6_r_EH;
B0 = k2_r_D - k2_r_E + R_25 * k5_r_DL - R_26 * k6_r_EH;

T_25_phi3 = [-sin(phi_25_3),-cos(phi_25_3),0;cos(phi_25_3),-sin(phi_25_3),0;0,0,0];
T_26_phi3 = [-sin(phi_26_3),-cos(phi_26_3),0;cos(phi_26_3),-sin(phi_26_3),0;0,0,0];

A = 2 * (T_25_phi3 * k5_r_DF)' * A0;
B = 2 * (T_26_phi3 * k6_r_EH)' * A0;
C = 2 * (T_25_phi3 * k5_r_DL)' * B0;
D = 2 * (T_26_phi3 * k6_r_EH)' * B0;

dphi_25_3_ds2 = 2 * s2 * B / (B * C - A * D);
dphi_26_3_ds2 = 2 * s2 * A / (B * C - A * D);

dphi_25_3_dt = dphi_25_3_ds2 * ds2;
dphi_26_3_dt = dphi_26_3_ds2 * ds2;

dR_25 = T_25_phi3 * dphi_25_3_dt;
dR_26 = T_26_phi3 * dphi_26_3_dt;

dT_25_phi3_dt = -R_25 * dphi_25_3_dt;
dT_26_phi3_dt = -R_26 * dphi_26_3_dt;

dA = 2 * (dT_25_phi3_dt * k5_r_DF)' * A0 ...
    + 2 * (T_25_phi3 * k5_r_DF)' * (dR_25 * k5_r_DF + dR_26 * k6_r_EH);
dB = 2 * (dT_26_phi3_dt * k6_r_EH)' * A0 ...
    + 2 * (T_26_phi3 * k6_r_EH)' * (dR_25 * k5_r_DF + dR_26 * k6_r_EH);
dC = 2 * (dT_25_phi3_dt * k5_r_DL)' * B0 ...
    + 2 * (T_25_phi3 * k5_r_DL)' * (dR_25 * k5_r_DL + dR_26 * k6_r_EH);
dD = 2 * (dT_26_phi3_dt * k6_r_EH)' * B0 ...
    + 2 * (T_26_phi3 * k6_r_EH)' * (dR_25 * k5_r_DL + dR_26 * k6_r_EH);

ddphi_25_3_ds2dt = 2 * ((ds2 * B + s2 * dB) * (B * C - A * D) ...
    - s2 * B * (dB * C + B * dC - dA * D - A * dD)) ...
    / ((B * C - A * D) ^ 2);
ddphi_26_3_ds2dt = 2 * ((ds2 * A + s2 * dA) * (B * C - A * D) ...
    - s2 * A * (dB * C + B * dC - dA * D - A * dD)) ...
    / ((B * C - A * D) ^ 2);

phi_5_3.phi_25_3 = phi_25_3;
phi_6_3.phi_26_3 = phi_26_3;
phi_5_3.dphi_5_3_ds2 = dphi_25_3_ds2;
phi_6_3.dphi_6_3_ds2 = dphi_26_3_ds2;
phi_5_3.ddphi_5_3_ds2dt = ddphi_25_3_ds2dt;
phi_6_3.ddphi_6_3_ds2dt = ddphi_26_3_ds2dt;
end

function tolerance = get_phi_3_func(phi_3,s2,L)
phi_25_3 = phi_3(1);
phi_26_3 = phi_3(2);

k2_r_D = L.k2_r_D;
k2_r_E = L.k2_r_E;
k5_r_DF = L.k5_r_F;
k5_r_DL = L.k5_r_L;
k6_r_EH = L.k6_r_H;

R_25 = [cos(phi_25_3),-sin(phi_25_3),0;sin(phi_25_3),cos(phi_25_3),0;0,0,1];
R_26 = [cos(phi_26_3),-sin(phi_26_3),0;sin(phi_26_3),cos(phi_26_3),0;0,0,1];

A0 = k2_r_D - k2_r_E + R_25 * k5_r_DF - R_26 * k6_r_EH;
B0 = k2_r_D - k2_r_E + R_25 * k5_r_DL - R_26 * k6_r_EH;

L_7 = L.L_7;
L_8 = L.L_8;
L_9 = L.L_9;

L_s2 = L_8 + L_9 + s2;
tolerance = (norm(A0) - L_7) ^ 2 + (norm(B0) - L_s2) ^ 2;
end