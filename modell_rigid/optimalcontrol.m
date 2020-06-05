

import casadi.*

T = 10; % Time horizon
N = 100; % number of control intervals

% Declare model variables
x = SX.sym('x',6);

x = [x(1); x(2); x(3); x(4); x(5); x(6)];
q= [x(1); x(2); x(3)]
dq= [x(4); x(5); x(6)]
u = SX.sym('u',3);
u = [u(1);u(2);u(3)];


% Model equations
xdot = [x(4);x(5);x(6);u(1);u(2);u(3)];                                    % 微分方程

% Objective term
L = (x(1)-pi/2)^2+(x(2)-2.6)^2+(x(3)-2.2)^2+5*u(1)^2+10*u(2)^2 +10*u(3)^2  ;        %评价函数

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
    lbw = [lbw; [-0.4;-0.4;-0.4]];              %对于u的约束？
    ubw = [ubw;  [0.3; 0.3; 0.3]];
    w0 = [w0;  0;0;0];

    % Integrate till the end of the interval
    Fk = F('x0',Xk,'p', Uk);
    Xk = Fk.xf;
    J=J+Fk.qf;

    % Add inequality constraint     %添加对于x的约束,
   
    g = {g{:}, Xk(1),Xk(2),Xk(3),Xk(4),Xk(5),Xk(6)};        %这里只约束了x1
    lbg = [lbg;  0; 0; 0;  0;  0;  0];
    ubg = [ubg;  pi/2; 2.6; 2.2;0.5;0.4;0.4];
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


% for  i = 1:100
  % q = x_opt(1:3,i)
  % dq = x_opt(4:6,i)
  %  [M,F] = Mass_Force_System(q,dq);
  %  pause(0.2);
 %   i=i+1
 %   F
%
%end

