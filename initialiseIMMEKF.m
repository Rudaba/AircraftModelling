function [X_IMM, P_IMM, Q, R_Noise, transMatrix, modeProbs] = initialiseIMMEKF(accelerations,omega)

%X_IMM and P_IMM represent the individual filter output
X_IMM           = zeros(30,6);
P_IMM           = zeros(30,30,6);
Q               = zeros(30,30,6);
R_Noise         = eye(6); %There'll only ever be one measurement

% X_IMM{1}        = zeros(30,1); 
% X_IMM{2}        = zeros(30,1);
% X_IMM{3}        = zeros(30,1);
% X_IMM{4}        = zeros(30,1);
% X_IMM{5}        = zeros(30,1);
% X_IMM{6}        = zeros(30,1);
% 
% P_IMM{1}        = zeros(30,30); 
% P_IMM{2}        = zeros(30,30);
% P_IMM{3}        = zeros(30,30);
% P_IMM{4}        = zeros(30,30);
% P_IMM{5}        = zeros(30,30);
% P_IMM{6}        = zeros(30,30);
% 
% Q{1}            = zeros(30,30); 
% Q{2}            = zeros(30,30);
% Q{3}            = zeros(30,30);
% Q{4}            = zeros(30,30);
% Q{5}            = zeros(30,30);
% Q{6}            = zeros(30,30);

%*****Model 1*****
X_IMM(1:3,1)   = [accelerations(1); accelerations(2); accelerations(3)]; % Body Acceleration
X_IMM(4:6,1)   = [omega(1); omega(2); omega(3)]; % Body Angular Rates

P_IMM(1,1,1)   = 1e-5^2;
P_IMM(2,2,1)   = 1e-5^2;
P_IMM(3,3,1)   = 1e-5^2;

P_IMM(4,4,1)   = 1e-5^2;
P_IMM(5,5,1)   = 1e-5^2;
P_IMM(6,6,1)   = 1e-5^2;

X_IMM(7,1)   = 1;  P_IMM(7,7,1)   = 1e-3^2; % CX_dE1
X_IMM(8,1)   = 1;  P_IMM(8,8,1)   = 1e-3^2; % CX_dE2

Q(1,1,1)       = 1e-5^2;
Q(2,2,1)       = 1e-5^2;
Q(3,3,1)       = 1e-5^2;

Q(4,4,1)       = 1e-5^2;
Q(5,5,1)       = 1e-5^2;
Q(6,6,1)       = 1e-5^2;

Q(7,7,1)       = 1e-5^2;
Q(8,8,1)       = 1e-5^2;

%*****Model 2*****
X_IMM(1:3,2)   = X_IMM(1:3,1);
X_IMM(4:6,2)   = X_IMM(4:6,1);

P_IMM(1,1,2)   = 1e-5^2;
P_IMM(2,2,2)   = 1e-5^2;
P_IMM(3,3,2)   = 1e-5^2;

P_IMM(4,4,2)   = 1e-5^2;
P_IMM(5,5,2)   = 1e-5^2;
P_IMM(6,6,2)   = 1e-5^2;

X_IMM(9,2)   = 1;  P_IMM(9,9,2)   = 1e-3^2; % CY_dR1
X_IMM(10,2)  = 1;  P_IMM(10,10,2) = 1e-3^2; % CY_dR2
X_IMM(11,2)  = 1;  P_IMM(11,11,2) = 1e-3^2; % CY_dE1

Q(1,1,2)       = 1e-5^2;
Q(2,2,2)       = 1e-5^2;
Q(3,3,2)       = 1e-5^2;

Q(4,4,2)       = 1e-5^2;
Q(5,5,2)       = 1e-5^2;
Q(6,6,2)       = 1e-5^2;

Q(9,9,2)       = 1e-5^2;
Q(10,10,2)     = 1e-5^2;
Q(11,11,2)     = 1e-5^2;

%*****Model 3*****
X_IMM(1:3,3)   = X_IMM(1:3,1);
X_IMM(4:6,3)   = X_IMM(4:6,1);

P_IMM(1,1,3)   = 1e-5^2;
P_IMM(2,2,3)   = 1e-5^2;
P_IMM(3,3,3)   = 1e-5^2;

P_IMM(4,4,3)   = 1e-5^2;
P_IMM(5,5,3)   = 1e-5^2;
P_IMM(6,6,3)   = 1e-5^2;

X_IMM(12,3) = 1;  P_IMM(12,12,3) = 1e-3^2; % CZ_dE1
X_IMM(13,3) = 1;  P_IMM(13,13,3) = 1e-3^2; % CZ_dE2
X_IMM(14,3) = 1;  P_IMM(14,14,3) = 1e-3^2; % CZ_dA1

Q(1,1,3)       = 1e-5^2;
Q(2,2,3)       = 1e-5^2;
Q(3,3,3)       = 1e-5^2;

Q(4,4,3)       = 1e-5^2;
Q(5,5,3)       = 1e-5^2;
Q(6,6,3)       = 1e-5^2;

Q(12,12,3)       = 1e-5^2;
Q(13,13,3)       = 1e-5^2;
Q(14,14,3)       = 1e-5^2;

%*****Model 4*****
X_IMM(1:3,4)   = X_IMM(1:3,1);
X_IMM(4:6,4)   = X_IMM(4:6,1);

P_IMM(1,1,4)   = 1e-5^2;
P_IMM(2,2,4)   = 1e-5^2;
P_IMM(3,3,4)   = 1e-5^2;

P_IMM(4,4,4)   = 1e-5^2;
P_IMM(5,5,4)   = 1e-5^2;
P_IMM(6,6,4)   = 1e-5^2;

X_IMM(15,4) = 1;  P_IMM(15,15,4) = 1e-3^2; % Cl_dA1
X_IMM(16,4) = 1;  P_IMM(16,16,4) = 1e-3^2; % Cl_dA2
X_IMM(17,4) = 1;  P_IMM(17,17,4) = 1e-3^2; % Cl_dA3
X_IMM(18,4) = 1;  P_IMM(18,18,4) = 1e-3^2; % Cl_dR1
X_IMM(19,4) = 1;  P_IMM(19,19,4) = 1e-3^2; % Cl_dR2
X_IMM(20,4) = 1;  P_IMM(20,20,4) = 1e-3^2; % Cl_dE1

Q(1,1,4)       = 1e-5^2;
Q(2,2,4)       = 1e-5^2;
Q(3,3,4)       = 1e-5^2;

Q(4,4,4)       = 1e-5^2;
Q(5,5,4)       = 1e-5^2;
Q(6,6,4)       = 1e-5^2;

Q(15,15,4)       = 1e-5^2;
Q(16,16,4)       = 1e-5^2;
Q(17,17,4)       = 1e-5^2;
Q(18,18,4)       = 1e-5^2;
Q(19,19,4)       = 1e-5^2;
Q(20,20,4)       = 1e-5^2;

%*****Model 5*****
X_IMM(1:3,5)   = X_IMM(1:3,1);
X_IMM(4:6,5)   = X_IMM(4:6,1);

P_IMM(1,1,5)   = 1e-5^2;
P_IMM(2,2,5)   = 1e-5^2;
P_IMM(3,3,5)   = 1e-5^2;

P_IMM(4,4,5)   = 1e-5^2;
P_IMM(5,5,5)   = 1e-5^2;
P_IMM(6,6,5)   = 1e-5^2;

X_IMM(21,5)  = 1;  P_IMM(21,21,5) = 1e-3^2; % Cm_dE1
X_IMM(22,5)  = 1;  P_IMM(22,22,5) = 1e-3^2; % Cm_dE2
X_IMM(23,5)  = 1;  P_IMM(23,23,5) = 1e-3^2; % Cm_dE3
X_IMM(24,5)  = 1;  P_IMM(24,24,5) = 1e-3^2; % Cm_dA1

Q(1,1,5)       = 1e-5^2;
Q(2,2,5)       = 1e-5^2;
Q(3,3,5)       = 1e-5^2;

Q(4,4,5)       = 1e-5^2;
Q(5,5,5)       = 1e-5^2;
Q(6,6,5)       = 1e-5^2;

Q(21,21,5)       = 1e-5^2;
Q(22,22,5)       = 1e-5^2;
Q(23,23,5)       = 1e-5^2;
Q(24,24,5)       = 1e-5^2;

%*****Model 6*****
X_IMM(1:3,6)   = X_IMM(1:3,1);
X_IMM(4:6,6)   = X_IMM(4:6,1);

P_IMM(1,1,6)   = 1e-5^2;
P_IMM(2,2,6)   = 1e-5^2;
P_IMM(3,3,6)   = 1e-5^2;

P_IMM(4,4,6)   = 1e-5^2;
P_IMM(5,5,6)   = 1e-5^2;
P_IMM(6,6,6)   = 1e-5^2;

X_IMM(25,6)  = 1;  P_IMM(25,25,6) = 1e-3^2; % Cn_dA1
X_IMM(26,6)  = 1;  P_IMM(26,26,6) = 1e-3^2; % Cn_dA2
X_IMM(27,6)  = 1;  P_IMM(27,27,6) = 1e-3^2; % Cn_dR1
X_IMM(28,6)  = 1;  P_IMM(28,28,6) = 1e-3^2; % Cn_dR2
X_IMM(29,6)  = 1;  P_IMM(29,29,6) = 1e-3^2; % Cn_dE1
X_IMM(30,6)  = 1;  P_IMM(30,30,6) = 1e-3^2; % Cn_dE2

Q(1,1,6)       = 1e-5^2;
Q(2,2,6)       = 1e-5^2;
Q(3,3,6)       = 1e-5^2;

Q(4,4,6)       = 1e-5^2;
Q(5,5,6)       = 1e-5^2;
Q(6,6,6)       = 1e-5^2;

Q(24,24,6)     = 1e-5^2;
Q(25,25,6)     = 1e-5^2;
Q(26,26,6)     = 1e-5^2;
Q(27,27,6)     = 1e-5^2;
Q(28,28,6)     = 1e-5^2;
Q(29,29,6)     = 1e-5^2;
Q(30,30,6)     = 1e-5^2;

%***Measurement Noise***
R_Noise(1,1)    = R_Noise(1,1)*(1e-2)^2;
R_Noise(2,2)    = R_Noise(2,2)*(1e-2)^2;
R_Noise(3,3)    = R_Noise(3,3)*(1e-2)^2;

transMatrix     = eye(6,6);
modeProbs       = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6];
