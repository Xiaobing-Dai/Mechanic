phi0 = 0;
s10 = -0.6;
s20 = 0;
phi = phi0;
s1 = s10;
s2 = s20;
u = [0;0;0];
for s1 = -0.6:0.01:0
    q = [phi;s1;s2];
    dq = zeros(3,1);
    [M,F] = Mass_Force_System(q,dq,u);
    pause(0.01);
end
s1 = 0;
for s2 = 0:0.01:1
    q = [phi;s1;s2];
    dq = zeros(3,1);
    [M,F] = Mass_Force_System(q,dq,u);
    pause(0.01);
end
s2 = 1;
for s1 = 0:0.01:2
    q = [phi;s1;s2];
    dq = zeros(3,1);
    [M,F] = Mass_Force_System(q,dq,u);
    pause(0.01);
end
s1 = 2;
for s2 = 1:0.01:2.2
    q = [phi;s1;s2];
    dq = zeros(3,1);
    [M,F] = Mass_Force_System(q,dq,u);
    pause(0.01);
end
s2 = 2.2;
for phi = 0:pi/20:2*pi
    q = [phi;s1;s2];
    dq = zeros(3,1);
    [M,F] = Mass_Force_System(q,dq,u);
    pause(0.01);
end
