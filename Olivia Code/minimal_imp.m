clear;
clc;
close all

% choose values for center of mass, plate geometry, fluid density, etc. All
% units cm g s

lcm = 0.1; % plate CoM
l = 2.54; % plate length
h = 0.05; % plate height
g = 980; % gravity
p = 1; % fluid density (water here)
m = 1; % plate 2D mass
I = .5; % plate 2D moment of inertia
mp = p * l * h; % mass of comparable fluid (for buoyancy-corrected mass)


%[x;y;theta;vxp;vyp;omega]--initial value in ODE 
x_0 = [0;0;-pi/6;5;0;0];

timerange = 5; % time to run sim

option=odeset('RelTol',1e-9, 'AbsTol',1e-8); % ODE options
F = ode(ODEFcn=@(t,x) Model(x,lcm,I,l,m,mp,p,g), InitialValue=x_0);
dt = 0.001; % timestep
Ntimes = timerange/dt;
t = linspace(0,timerange,Ntimes);

sol = solve(F,t); % ODE solve

% extract solution
t = sol.Time;
solution = sol.Solution;

x = solution(1,:);
y = solution(2,:);

plot(x,y);
xlabel('x');
ylabel('y');









function y = Model(initial,lcm,I,l,m,mp,p,g) 
    % y(1:6)= dot_x, dot_y, dot_theta, dot_vxp, dot_vyp, dot_omega
    % x&y are the geometry center, x prime&y prime are center of mass
    y = zeros(size(initial)); 

    CL1 = 5.2;
    CL2 = 0.95;
    CD0 = 0.1;
    CD1 = 5;
    CD2 = 1.9;
    CR = 1.1;
    CCP0 = 0.3;
    CCP1 = 3.5;
    CCP2 = 0.2;

    
    m1 = 0;                      % added mass m_11
    m2 = pi * p * l^2 * 0.25;    % added mass m_22
    Ia = pi * p *l^4* (1+32*lcm^2/l^2)/128; % added moment of inertia
    
    c = cos(initial(3));             % cos(theta)
    s = sin(initial(3));             % sin(theta)
    vxp = initial(4);                % Vx prime
    vyp = initial(5);                % Vy prime
    omega = initial(6);                % omega, rotational velocity
    
    v_cg = sqrt(vxp^2+(vyp-omega*lcm)^2); % velocity of geometric center
    alpha = atan2(vyp-omega*lcm,vxp);  %angle of attack
    
    
    %weighting function
    as = 14*pi/180;
    d =  6*pi/180;
    f1 = ( 1 - tanh((abs(alpha)-as)/d) )/2;   % for |alpha| < pi/2
    f2 = ( 1 - tanh((pi-abs(alpha)-as)/d) )/2;% for |alpha| > pi/2
    
    %lift coeff
    CT_1= CL1*f1*sin(alpha)+ (1-f1)*CL2*sin(2*alpha); %(-90,90)
    CT_2= -(CL1*f2*(pi-abs(alpha))+(1-f2)*CL2*sin(2*(pi-abs(alpha))));%(90,180)
    CT_3= CL1*f2*(pi-abs(alpha))+ (1-f2)*CL2*sin(2*(pi-abs(alpha)));%(-180,-90)
    
    %drag coeff
    CD_1 = f1*(CD0 + CD1*(sin(alpha))^2) + (1-f1)*CD2*(sin(alpha))^2;
    CD_2=f2*(CD0 + CD1*(pi-abs(alpha))^2) + (1-f2)*CD2*(sin(pi-alpha))^2;
    
    %center of pressure
    r=abs(alpha);
    rr=pi-r;
    lcp_1=f1*(CCP0-CCP1*r^2)+(1-f1)*(CCP2 -(CCP2/(pi/2))*r);
    lcp_1=lcp_1*l;
    lcp_2=-(f2*(CCP0-CCP1*rr^2)+(1-f2)*(CCP2 -(CCP2/(pi/2))*rr));
    lcp_2=lcp_2*l;
    
    if r<=pi/2
    CT_alpha=CT_1;
    CD_alpha=CD_1;
    l_cp=lcp_1;    
    elseif alpha>pi/2
    CT_alpha=CT_2;
    CD_alpha=CD_2;
    l_cp=lcp_2; 
    else
    CT_alpha=CT_3;
    CD_alpha=CD_2;
    l_cp=lcp_2;
    end
    
    l_tau = l_cp - lcm;
    
    lcr = 0;
    
    %lift force   
    L_0 =0.5*p*l*(CT_alpha*v_cg-l*CR*omega);    
    Lx=L_0*(vyp-omega*lcm);
    Ly=-L_0*vxp;
    
    % drag force
    D_0 = -0.5*p*l*CD_alpha*v_cg; 
    D_xp=D_0*vxp;
    D_yp=D_0*(vyp-omega*lcm);
    
    
    % torques
    Tt = -0.5*p*l*v_cg*(CT_alpha*vxp+CD_alpha*(vyp-omega*lcm))*l_tau; % torque from translational lift
    
    Tcr = -0.5*p*l^2*CR*omega*vxp*(lcm-lcr); % rotational torque
    Tr= -1/128*p*l^4*CD2*omega*abs(omega)*((2*lcm/l+1)^4+(2*lcm/l-1)^4);%damping torque
    Tb=-mp*g*lcm*c; % buoyancy
    
    %ODE function
    total_torque= Tt+Tr+Tb + Tcr;   %total torque
    F_xp = (m+m2)*omega*vyp -m2*omega^2*lcm+ Lx + D_xp - (m-mp)*g*s ;  % total forces in xp-"F_xp"
    F_yp = -(m+m1)*omega*vxp+m2*total_torque/(I+Ia)*lcm+ Ly + D_yp - (m-mp)*g*c ;% total forces in yp_"F_yp"
    
    
    y(1) = vxp*c-vyp*s;
    y(2) = vxp*s+vyp*c;
    y(3) = omega;
    y(4) = F_xp / (m+m1);
    y(5) = F_yp / (m+m2);
    y(6) = total_torque / (I+Ia);  
end