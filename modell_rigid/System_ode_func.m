function dx = System_ode_func(t,x)
q = x(1:3);
dq = x(4:6);

[M,F] = Mass_Force_System(q,dq);
ddq = M\F;

dx = [dq;ddq];

end