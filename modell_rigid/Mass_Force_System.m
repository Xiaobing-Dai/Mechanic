function [M,F] = Mass_Force_System(q,dq,u)
% clc;
% clear;
% phi = 0;
% s1 = 2;
% s2 = 0;
% q = [phi;ds1;ds2];
% s1
% dq = zeros(3,1);
% dq(1) = 5;
% dq(2) = 4;
% dq(3) = 3;

% m_tot = 100;
% rho = 1;
% k_theta_0 = eye(3);
% theta(1,1) = 100;
% theta(2,2) = 33.3;
% theta(3,3) = 133.33;
% k_r_0c = [0.666666666666667,0.866,0]';
% u = zeros(6,1);
% [M_e,F_e,others] = Mass_Force_element(q,dq,m_tot,rho,k_theta_0,k_r_0c,u)

%%
L.k1_r_0A = [-0.866,1.5,0]';
L.k1_r_0B = [0.866,0.5,0]';
L.L_A0A = 5;
L.L_3 = 2;
L.k3_AB = [-1.5,0,0]';
L.L_4 = 2;
L.k3_B0BB = [1.5,0,0]';
L.k2_r_A = [5,0,0]';
L.k2_r_D = [13,0,0]';
L.k2_r_E = [12.5,-1,0]';

L.k5_r_F = [0.5,-1,0]';
L.k5_r_L = [5.5,0,0]'; %
L.k5_r_P = [13,0,0]';
L.L_6 = 1.5;
L.k6_r_H = [1.5,0,0]';
L.L_7 = 2;
L.L_8 = 2.18;
L.k8_HGG = [1.5,0,0]';  %
L.L_9 = 2.18;
L.k9_LG = [-1.5,0,0]'; %
L.k7_r_H = [2,0,0]';
%%
[M_1,F_1,Body{1},Point_0] = Body1(q,dq,u,L);
Point_A0 = Point_on_Body(Body{1},Point_0,L.k1_r_0A);
Point_B0 = Point_on_Body(Body{1},Point_0,L.k1_r_0B);
[M_2,F_2,Body{2},phi_2_3] = Body2(q,dq,u,{Body{1},Point_A0,Point_B0},L);
Point_A = Point_on_Body(Body{2},Point_A0,L.k2_r_A);
[M_4,F_4,Body{4}] = Body4(q,dq,u,{Body{1},Point_A0,Point_B0},L);
[M_3,F_3,Body{3}] = Body3(q,dq,u,{Body{2},Body{4},Point_A},L);
% M_1 + M_2 + M_3 + M_4
%%
Point_B = Point_on_Body(Body{3},Point_A,L.k3_AB);
Point_BB = Point_on_Body(Body{4},Point_B0,L.k3_B0BB);
%%
Point_D = Point_on_Body(Body{2},Point_A0,L.k2_r_D);
Point_E = Point_on_Body(Body{2},Point_A0,L.k2_r_E);
[phi_5_3,phi_6_3] = get_phi_up_56(q,dq,L);
[M_5,F_5,Body{5}] = Body5(q,dq,u,{Body{2},Point_D,phi_5_3},L);
Point_F = Point_on_Body(Body{5},Point_D,L.k5_r_F);
Point_L = Point_on_Body(Body{5},Point_D,L.k5_r_L);
Point_P = Point_on_Body(Body{5},Point_D,L.k5_r_P);
[M_6,F_6,Body{6}] = Body6(q,dq,u,{Body{2},Point_E,phi_6_3},L);
Point_H = Point_on_Body(Body{6},Point_E,L.k6_r_H);
[phi_7_3,phi_8_3,phi_9_3] = get_phi_up_789(q,dq,L,{Body{2},Body{6},Point_D,Point_H,phi_2_3});
[M_7,F_7,Body{7}] = Body7(q,dq,u,{Body{6},Point_H,phi_7_3},L);
[M_8,F_8,Body{8}] = Body8(q,dq,u,{Body{6},Point_H,phi_8_3},L);
[M_9,F_9,Body{9}] = Body9(q,dq,u,{Body{5},Point_L,phi_9_3},L);
%%
Point_G = Point_on_Body(Body{9},Point_L,L.k9_LG);
Point_GG = Point_on_Body(Body{8},Point_H,L.k8_HGG);
%%
M = M_1 + M_2 + M_3 + M_4 + M_5 + M_6 + M_7 + M_8 + M_9;
F = F_1 + F_2 + F_3 + F_4 + F_5 + F_6 + F_7 + F_8 + F_9;
%% ODE Mddq-F=0
Print_Body1 = [zeros(3,1),Point_A0.r,Point_B0.r,zeros(3,1)];
Print_Body2 = [Point_A0.r,Point_D.r,Point_E.r,Point_A.r];
% Print_Body34 = [Point_A.r,Point_B0.r];
Print_Body3 = [Point_A.r,Point_B.r];
Print_Body4 = [Point_B0.r,Point_BB.r];
Print_Body5 = [Point_D.r,Point_F.r,Point_L.r,Point_P.r,Point_D.r];
Print_Body6 = [Point_E.r,Point_H.r];
Print_Body7 = [Point_F.r,Point_H.r];
% Print_Body89 = [Point_L.r,Point_H.r];
Print_Body8 = [Point_H.r,Point_GG.r];
Print_Body9 = [Point_L.r,Point_G.r];
hold off;
plot3(0,0,0);

hold on;
plot3(Print_Body1(1,:),Print_Body1(2,:),Print_Body1(3,:),'k-');

h=plot3(Print_Body2(1,:),Print_Body2(2,:),Print_Body2(3,:),'b-');

% plot3(Print_Body34(1,:),Print_Body34(2,:),Print_Body34(3,:),'g-');
plot3(Print_Body3(1,:),Print_Body3(2,:),Print_Body3(3,:),'g-');
plot3(Print_Body4(1,:),Print_Body4(2,:),Print_Body4(3,:),'g-');
plot3(Print_Body5(1,:),Print_Body5(2,:),Print_Body5(3,:),'b-');
plot3(Print_Body6(1,:),Print_Body6(2,:),Print_Body6(3,:),'r-');
plot3(Print_Body7(1,:),Print_Body7(2,:),Print_Body7(3,:),'r-');
% plot3(Print_Body89(1,:),Print_Body89(2,:),Print_Body89(3,:),'g-');
plot3(Print_Body8(1,:),Print_Body8(2,:),Print_Body8(3,:),'g-');
plot3(Print_Body9(1,:),Print_Body9(2,:),Print_Body9(3,:),'g-');
plot3([1.7322,-3,-3],[0,2.7319,0],[0,0,0]);
    


grid on
 view(180,300)
set(gca,'xdir','reverse')
axis equal;
xlim([-20,20]);
ylim([0,30]);
zlim([-10,10]);



