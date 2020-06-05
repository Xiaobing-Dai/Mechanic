clear;
tspan = [0 5];
x0 = [0,2,2,0,0,0]';
[t,x] = ode45(@(t,x)Rigid_Body_odefun_Test(t,x),[0 10],x0)

function [dxdt] = Rigid_Body_odefun_Test(t,x)
q = x(1:3);
dq = x(4:6);
[M,F] = Mass_Force_System(q,dq);

ddq = M\F
dxdt = [dq;ddq];
pause(0.01);
end