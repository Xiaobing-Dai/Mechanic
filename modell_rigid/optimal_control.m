function [dxxdt] = optimal_control(t,xx)
% xx = q dq phi = q u phi
% dxxdt = dq ddq dphi = u ddq dphi
% Guetemass
clc;
% % clear all;
J = 1;  %改成正常的
% % Dimension von G?tema? ist immer 1x1. Hier G?tema? eingeben.
y = [1;1;1;1;1;1];
u = [0;0;0];  % input [3x1]
x = [6;5;4;3;2;1]; %Zustand [6x1]

dxxdt = zeros(numel(xx),1);% dq ddq dphy

x = xx(1:6);
q = x(1:3);
dq = x(4:6);

% x=[q,dq];%q 3x1 x 6x1

[M,F] = Mass_Force_System(q,dq);

ddq=M\F;
dxdt = [dq;ddq];

% M = get_M(x,u);  % 应该直接从动力学方程中倒入
% F = get_F(x,u);
% M_F = M\F;  %[3x1]
g = get_g(u,x,y);     % g in form function
f = [dq ; ddq];
% phy = [6;5;4;3;2;1];
phy = xx(7:12);

% H = g + phy'*f ;      %phy ist lagelangri suan zi, f ist Dynamik

%minimum H in Abh?ngigkeit von u, weil u begrenzt ist
u0 = [0;0;0]; 
[u,~] = fmincon(@(uu)get_H(uu,x,y,phy),u0,[],[],[],[],[0,0,0],[0.3,0.3,0.3],[]);
% dxdt = dq ddq = u ddq
dxdt(1:3) = u;

%%
dgdx = zeros(numel(x),1); %l?sen dgdx [6x1] 
h = 0.01;
for i = 1 : numel(dgdx)
    dgdx(i)=get_dgdx(u,x,y,i,h); 
end
%%
phyf = phy(1)*f(1) + phy(2)*f(2)+ phy(3)*f(3)+...
    phy(4)*f(4)+ phy(5)*f(5)+ phy(6)*f(6);
dphyfdx = zeros(numel(x),1);
h = 0.01;
for i = 1 : numel(dphyfdx)
    dphyfdx(i) = get_dphydx(u,x,y,i,h,phy,f);
end

dHdx = dgdx + dphyfdx;
dphydt = -dHdx;      %H de pian dao

%minimum H in Abh?ngigkeit von u, weil u begrenzt ist
% [uu,H] = fmincon(@(uu)get_H(uu,x,y,phy),u,[],[],[],[],[-100,-100,-100],[100,100,100],[]);  

% 每个时刻时刻都是取让H最小的u 
% u soll hier begrenzt werden

% xx = [x,phy]; % phi 6x1 xx 12x1
dxxdt = [dxdt;dphydt];
end
% tspan = [0 5];
% x0 = [0,-0.6,-0.7,0,0,0,0,0,0,0,0,0]
% [t,xx] = ode45(optimal_control,tspan,x0)



function H = get_H(u,x,y,phy)
% M = get_M(x,u);
% F = get_F(x,u);
 q = x(1:3);
%  dq = x(4:6);
 dq = u;
 [M,F] = Mass_Force_System(q,dq);
 ddq=M\F;
 f = [u;ddq];
% M_F = M\F;
% f = [u ; M_F];
g = get_g(u,x,y);
H = g + phy'*f ;
end


% function M = get_M(x,uu)
% M = magic(3);
% end
% 
% function F = get_F(x,uu)
% F = ones(3,1);
% end

function dgdx = get_dgdx(u,x,y,n,h)
h_set = [-2,-1,1,2] * h;
g_set = zeros(1,numel(h_set));

for i = 1:numel(h_set)
    delta_h = h_set(i);
    x(n) = x(n) + delta_h;  
    g_set(i) = get_g(u,x,y);
end

dgdx = 1/(12*h) * (g_set(1) - 8 * g_set(2) + 8 * g_set(3) - g_set(4));
end

function g=get_g(u,x,y)
xe = [0;2;2;0;0;0];

deltax= x-xe;
g = 1/2 * deltax' * eye(6)* deltax;  %不'同的评价函数主要修改g

end


function dphyfdx = get_dphydx(u,x,y,n,h,phy,f)
h_set = [-2,-1,1,2] * h;
g_set = zeros(1,numel(h_set));

for i = 1:numel(h_set)
    delta_h = h_set(i);
    x(n) = x(n) + delta_h;
%     f(x) = [dq,ddq]
    q = x(1:3);
    dq = x(4:6);
    [M,F] = Mass_Force_System(q,dq);
    ddq = M\F;
    f = [dq;ddq];
    g_set(i) = get_phyf(phy,f);
end

dphyfdx = 1/(12*h) * (g_set(1) - 8 * g_set(2) + 8 * g_set(3) - g_set(4));
end

function phyf = get_phyf(phy,f)
phyf = phy(1)*f(1) + phy(2)*f(2)+ phy(3)*f(3)+...
    phy(4)*f(4)+ phy(5)*f(5)+ phy(6)*f(6);
end

import casadi.*

T = 10; % Time horizon
N = 200; % number of control intervals

% Declare model variables
x = SX.sym('x',6);

x = [x(1); x(2); x(3); x(4); x(5); x(6)];
u = SX.sym('u',3);
u = [u(1);u(2);u(3)];


% Model equations
xdot = [x(4);x(5);x(6);u(1);u(2);u(3)];                                    % 微分方程

% Objective term
L = (x(1)-2)^2+ (x(2)-1)^2+ (x(3)-3)^2+ (x(4)-0.4)^2+ (x(5)-0.5)^2+ (x(6)-0.6)^2 ;        %评价函数

% Formulate discrete time dynamics
if false
   % CVODES from the SUNDIALS suite
   dae = struct('x',x,'p',u,'ode',xdot,'quad',L);    %求ode
   opts = struct('tf',T/N);
   F = integrator('F', 'cvodes', dae, opts);
else
   % Fixed step Runge-Kutta 4 integrator
   M = 4; % RK4 steps per interval
   DT = T/N/M;
   f = Function('f', {x, u}, {xdot, L});
   X0 = MX.sym('X0',6);     %注意维度
   U = MX.sym('U',3);         %注意维度
   X = X0;
   Q = 0;
   for j=1:M
       [k1, k1_q] = f(X, U);
       [k2, k2_q] = f(X + DT/2 * k1, U);
       [k3, k3_q] = f(X + DT/2 * k2, U);
       [k4, k4_q] = f(X + DT * k3, U);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
    end
    F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'}); %对应95行
end

% Evaluate at a test point
Fk = F('x0',[0.2; 0.3;0.3;0.3;0.3;0.3],'p',0.4);
disp(Fk.xf)
disp(Fk.qf)

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% Formulate the NLP
Xk = [0;0;0;0;0;0];                        %Xk是什么: Initial conditions
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)],1,3);
    w = {w{:}, Uk(:)};
    lbw = [lbw; [-0.5;-0.5;-0.5]];              %对于u的约束？
    ubw = [ubw;  [0.5; 0.5; 0.5]];
    w0 = [w0;  0;0;0];

    % Integrate till the end of the interval
    Fk = F('x0',Xk,'p', Uk);
    Xk = Fk.xf;
    J=J+Fk.qf;

    % Add inequality constraint     %添加对于x的约束,
   
    g = {g{:}, Xk(1),Xk(2),Xk(3),Xk(4),Xk(5),Xk(6)};        %这里只约束了x1
    lbg = [lbg;  0; 0; 0;  0;  0;  0];
    ubg = [ubg;  2; 3; 4;0.4;0.4;0.4];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
             'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

% Plot the solution
u_opt = w_opt;
x_opt = [0;0;0;0;0;0];                  % 画出 x1,x2的初始点       
u = [];
for k=0:N-1
    Fk = F('x0', x_opt(:,end), 'p', u_opt(3*k+1:3*k+3));
    x_opt = [x_opt, full(Fk.xf)];
    u = [u,[u_opt(3*k+1:3*k+3)]];
end
u = u';

x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
x3_opt = x_opt(3,:);
x4_opt = x_opt(4,:);
x5_opt = x_opt(5,:);
x6_opt = x_opt(6,:);
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x1_opt, '--')
plot(tgrid, x2_opt, '--')
plot(tgrid, x3_opt, '--')
plot(tgrid, x4_opt, '-')
plot(tgrid, x5_opt, '-')
plot(tgrid, x6_opt, '-')
stairs(tgrid, [u(:,1); nan], '-.')
stairs(tgrid, [u(:,2); nan], '-.')
stairs(tgrid, [u(:,3); nan], '-.')
xlabel('t')
legend('x1','x2','x3','x4','x5','x6','u1','u2','u3')

