syms dE dA dR throttle 
syms qbar ALP BET Vt dt H
syms ax ay az p q r  
syms prev_p prev_q prev_r
syms CX_dE 
syms CY_dR CY_dE
syms CZ_dE CZ_dA
syms Cl_dA Cl_dR Cl_dE
syms Cm_dE Cm_dA
syms Cn_dA Cn_dR Cn_dE

rho     = 1.225;
vs      = 340.3;
Ix      = 24970*14.593903*0.092903/15;
Iy      = 122190*14.593903*0.092903/15;
Iz      = 139800*14.593903*0.092903/15;
Ixz     = 1175*14.593903*0.092903/15;
b       = 20/3;
cbar    = 3;
S       = 20;
m       = 38924 * 0.453592 / 15;

xCG_ref = 0;
xCG     = 0;

hT      = H/3048;

Tmax = ((30.21-0.668*hT-6.877*hT^2+1.951*hT^3-0.1512*hT^4) + ...
    (Vt/vs).*(-33.8+3.347*hT+18.13*hT^2-5.865*hT^3+0.4757*hT^4) + ...
    (Vt/vs)^2.*(100.8-77.56*hT+5.441*hT^2+2.864*hT^3-0.3355*hT^4) + ...
    (Vt/vs)^3.*(-78.99+101.4*hT-30.28*hT^2+3.236*hT^3-0.1089*hT^4) + ...
    (Vt/vs)^4.*(18.74-31.6*hT+12.04*hT^2-1.785*hT^3+0.09417*hT^4))*4448.22/20; % Newton's

T = Tmax*throttle;

CX_dE1 = CX_dE*9.5e-4; 
CX_dE2 = CX_dE*8.5e-7; 

CY_dR1 = CY_dR*1.55e-3;
CY_dR2 = CY_dR*8e-6;
CY_dE1 = CY_dE*1.75e-4;

CZ_dE1 = CZ_dE*4.76e-3; 
CZ_dE2 = CZ_dE*3.3e-5; 
CZ_dA1 = CZ_dA*7.5e-5;

Cl_dA1 = Cl_dA*6.1e-4;
Cl_dA2 = Cl_dA*2.5e-5;
Cl_dA3 = Cl_dA*2.6e-6;
Cl_dR1 = Cl_dR*-2.3e-4;
Cl_dR2 = Cl_dR*4.5e-6;
Cl_dE1 = Cl_dE*5.24e-5; 

Cm_dE1 = Cm_dE*6.54e-3; 
Cm_dE2 = Cm_dE*8.49e-5; 
Cm_dE3 = Cm_dE*3.74e-6;
Cm_dA1 = Cm_dA*3.5e-5;

Cn_dA1 = Cn_dA*1.4e-5;
Cn_dA2 = Cn_dA*7.0e-6;
Cn_dR1 = Cn_dR*9.0e-4; 
Cn_dR2 = Cn_dR*4.0e-6; 
Cn_dE1 = Cn_dE*8.73e-5;
Cn_dE2 = Cn_dE*8.7e-6;

CX = -0.0434 + 2.39e-3*ALP+2.53e-5*BET^2-1.07e-6*ALP*BET^2+CX_dE1*dE-CX_dE2*dE*BET^2+(180*q*cbar/pi/2/Vt)*(8.73e-3+0.001*ALP-1.75e-4*ALP^2);
CY = -0.012*BET+CY_dR1*dR-CY_dR2*dR*ALP+(180*b/pi/2/Vt)*(2.25e-3*p+0.0117*r-3.67e-4*r*ALP+CY_dE1*r.*dE);
CZ = -0.131-0.0538*ALP-CZ_dE1*dE-CZ_dE2*dE*ALP-CZ_dA1*dA.^2+(180*q*cbar/pi/2/Vt).*(-0.111+5.17e-3*ALP-1.1e-3*ALP^2);
Cl = -5.98e-4*BET-2.83e-4*ALP*BET+1.51e-5*ALP^2*BET-dA.*(Cl_dA1+Cl_dA2*ALP-Cl_dA3*ALP^2)-dR.*(Cl_dR1+Cl_dR2*ALP)+(180*b/2/pi/Vt).*(-4.12e-3*p-5.24e-4*p*ALP+4.36e-5*p*ALP^2+4.36e-4*r+1.05e-4*r*ALP+Cl_dE1*r.*dE);
Cm = -6.61e-3-2.67e-3*ALP-6.48e-5*BET^2-2.65e-6*ALP*BET^2-Cm_dE1*dE-Cm_dE2*dE*ALP+Cm_dE3*dE*BET^2-Cm_dA1*dA.^2+(180*q*cbar/2/pi/Vt).*(-0.0473-1.57e-3*ALP)+(xCG_ref-xCG)*CZ;
Cn = 2.28e-3*BET+1.79e-6*BET^3+Cn_dA1*dA+Cn_dA2*dA*ALP-Cn_dR1*dR+Cn_dR2*dR*ALP+(180*b/2/pi/Vt).*(-6.63e-5*p-1.92e-5*p*ALP+5.06e-6*p*ALP^2-6.06e-3*r-Cn_dE1*r.*dE+Cn_dE2*r.*dE*ALP)-cbar/b*(xCG_ref-xCG)*CY;

c1 = ((Iy-Iz)*Iz-Ixz^2)/(Ix*Iz-Ixz^2);
c2 = (Ix-Iy+Iz)*Ixz/(Ix*Iz-Ixz^2);
c3 = Iz/(Ix*Iz-Ixz^2);
c4 = Ixz/(Ix*Iz-Ixz^2);
c5 = (Iz-Ix)/Iy;
c6 = Ixz/Iy;
c7 = 1/Iy;
c8 = ((Ix-Iy)*Ix-Ixz^2)/(Ix*Iz-Ixz^2);
c9 = Ix/(Ix*Iz-Ixz^2);

heng = 0;
pdot = (c1*r+c2*p+c4*heng)*q+qbar*S*b*(c3*Cl+c4*Cn);
qdot = (c5*p-c7*heng)*r-c6*(p^2-r^2)+qbar*S*cbar*c7*Cm;
rdot = (c8*p-c2*r+c9*heng)*q+qbar*S*b*(c4*Cl+c9*Cn);

%Equations
axNew = (qbar*S*CX+T)/m; 
ayNew = qbar*S*CY/m; 
azNew = qbar*S*CZ/m;

pNew = p + dt*(prev_p+pdot)/2;
qNew = q + dt*(prev_q+qdot)/2;
rNew = r + dt*(prev_r+rdot)/2;

stateVec = [axNew; ayNew; azNew;...
    pNew;qNew;rNew;...
    CX_dE;...
CY_dR; CY_dE;...
CZ_dE; CZ_dA;...
Cl_dA; Cl_dR; Cl_dE;...
Cm_dE; Cm_dA;...
Cn_dA; Cn_dR; Cn_dE];

states   = [ax; ay; az;...
    p;q;r;...
    CX_dE;...
CY_dR; CY_dE;...
CZ_dE; CZ_dA;...
Cl_dA; Cl_dR; Cl_dE;...
Cm_dE; Cm_dA;...
Cn_dA; Cn_dR; Cn_dE];

F = vpa(jacobian(stateVec,states));

