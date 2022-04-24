% Implementation of Extended kalman filter according ot the paper:
% "Extended Range PMSM Sensorless Speed Drive Based on Stochastic Filtering"

% by Nicola Chiodetto @08/06/2020

close all; clear; clc;

%% MACHINE AND ELECTRONICS PARAMETERS

% parameters from Bolognani paper of EKF electro-magnetic
Ld      = 3e-3;         % d-axis inductance 
Lq      = 3e-3;         % q-axis inductance
R       = 1.9;          % phase resistance (this might change real time due to transient)
phi_PM  = 0.1;          % PM flux (this might change during transient as well)
p       = 4;            % pole pairs
rpm     = 4000;         % rated speed
ff      = rpm*p/60;     % fundamental frequency
wm      = 2*pi*ff;      % ang elec speed
fpwm    = 20e3;         % switching freq Hz
Ts      = 1/fpwm;       % switching period (for now assumed equal to sampling)
Vdc     = 450;          % DC-Link voltage

% mechanical
inertia     = 1.8e-4;
LT          = 2.8;      % load torque
LT_time     = 0.15;     % load torque kick-in time
ss_speed    = 4/5;      % scaling start speed (to reduce computation transient)

% Fault parameters and power electronics parameters (other stuff )
Load_R  = 0.1;      %Short circuit test
OC_R    = 100000;   %Open  circuit test
Cbus    = 100e-6;   % bus capacitance
Rbus    = 0.001;    % bus cap ESR

%% Speed sensor features - AD2S1210 

% All the speed sensor features are grouped into structure ssensor

ssensor.resolution  = 12;           % Resolution choice
ssensor.t           = 1/4.096e6;    % sampling period
ssensor.k1          = 1.8/2.5;      % ADC gain; valid for all cases

switch ssensor.resolution 
    
    case 10    
        ssensor.k2  = 6*10^6 * 2*pi; % error gain
        ssensor.a   = 8187/8192; % compensator zero coeff
        ssensor.b   = 509/512; % compensator pole coeff
        ssensor.c   = 1/1024000; % integrator gain
        
    case 12        
        ssensor.k2  = 18*10^6 * 2*pi; % error gain
        ssensor.a   = 4095/4096; % compensator zero coeff
        ssensor.b   = 4085/4096; % compensator pole coeff
        ssensor.c   = 1/4096000; % integrator gain
        
    case 14        
        ssensor.k2  = 82*10^6 * 2*pi; % error gain
        ssensor.a   = 8191/8192; % compensator zero coeff
        ssensor.b   = 16359/16384; % compensator pole coeff
        ssensor.c   = 1/16384000; % integrator gain
        
    case 16        
        ssensor.k2  = 66*10^6 * 2*pi; % error gain
        ssensor.a   = 32767/32768; % compensator zero coeff
        ssensor.b   = 32757/32768; % compensator pole coeff
        ssensor.c   = 1/65536000; % integrator gain
        
    otherwise
        error('Resolution not permitted.')
        
end

% parameters reduced order transfer function ------------------------------
t1 = ssensor.t*(1+ssensor.a)/2/(1-ssensor.a);
t2 = ssensor.t*(1+ssensor.b)/2/(1-ssensor.b);
Ka = ssensor.k1*ssensor.k2*(1-ssensor.a)/(ssensor.a-ssensor.b);

% Complete transfer function ----------------------------------------------
% numerator coeff
a3 = ssensor.t^3*( 1 + ssensor.a );
a2 = ssensor.t^2*( -2*ssensor.a + 2 + 4 + 4*ssensor.a );
a1 = ssensor.t*( 4+4*ssensor.a-8*ssensor.a+8 );
a0 = -8*ssensor.a + 8;

% denominator coeff
b3 = ssensor.t^3*(1+ssensor.b);
b2 = ssensor.t^2*(2-2*ssensor.b);

K = ssensor.c^2*ssensor.k1*ssensor.k2/4; % gain

G = K * tf([a3 a2 a1 a0],[b3 b2 0 0]); % cont time tf - Open Loop
H = G / (1+G); % cont time tf - Closed Loop
[num den] = tfdata(H, 'v'); % extract num and den

% augment of 1 degree for numerator for the speed
num_speed = [num 0];

%% NOISE implementation
% assuming noise to be normal (generate the Gauss noise)
snoise_mean = 0;
snoise_std = 0.5;
switch ssensor.resolution
    case 12
        gainSpeed_noise = 29.3 /60*2*pi;
        gainAngle_noise = 0.09/180*pi;
    case 14
        gainSpeed_noise = 3.6 /60*2*pi;
        gainAngle_noise = 0;
end

%% CURRENT ERROR SENSOR MODELING
imean       = 0;
iaccuracy   = 12.8* (40+40)/100;
istdDev     = iaccuracy/3;
ivar        = istdDev^2; 

%% MODELLING 

% model for the 2DOF - IMC internal model control (blue box in current control block in SImulink model)
% general matrices
I = eye(2);         % indentity matrix
J = skewdec(2,0);   % Orthogonal rotation matrix: this is a skew-symmetric matrix

% Resistance Matrices
Rs  = R;%*eye(2);
b   = [Rs/Ld; 0];

% Inductance Matrices
L       = [Ld 0 ;0 Lq];
Linv    = inv(L); 

% Permanent magnet flux vector
d = [-1/Ld; 0];

% stator fluxes as state variables (4 state space feedback)
A = [-Rs/Ld wm;
    -wm -Rs/Lq];
C = [1/Ld 0;
    0 1/Lq];

%% CONTROLLERS

% internal model control is implemented with back calculation anti-windup
alpha = 2*pi*200;   % chosen bandwidth for close loop emulating low pass filter
% matrices useful for control
Bc  = alpha*Rs*eye(2);
Dc  = alpha*[Ld 0; 0 Lq];
W   = [0 -Lq; Ld 0];    % decoupling needs to be adjusted with speed multiplication
K   = Dc + Ts*Bc;       % proportional control matrix
Ki  = Ts*Bc*inv(K);     % integral control matrix

%% SPEED CONTROLLER -------------------------------------------------------
% The speed controller is based on:
% "A General Algorithm for Speed and Position Estimation of AC Motor"

Tspeed  = Ts*10/10; % by Tspeed we decide the period of the speed loop 
T_nom   = LT;
t_load  = 0.3; % to be put to zero for the transient start up from 0 to 6krpm

T_lim   = 5;    % saturation antiwindup for the torque reference
P_lim   = 1.5e3; % imposing a limit in the power and getting the torque out

t_ro        = inertia*wm/p/T_nom*2;
t_rc        = 0.5*t_ro;
alphas      = log(9)/t_rc;
zeta        = 0.9;

% gains
Kp_speed    = alphas*inertia; 
Ki_speed    = (alphas/2/zeta)^2*inertia;

Kbs         = 2/Kp_speed; % anti-windup back calculation tuning parameter

% plant and controller transfer functions
G = tf([1],[inertia 0]); 
F = tf([Kp_speed Ki_speed],[1 0]);

% standard 2DOF
Fr = tf([inertia*(alphas/2/zeta)^2],[1 alphas]); 

% servo-acceleration 2DOF
m       = 3;
J_coeff = inertia*alphas;
Fr      = tf([J_coeff*(2*m-1) J_coeff*alphas*(m^2-1/(4*zeta^2)) 0],[1 2*m*alphas (m*alphas)^2]); 

Lag     = tf([-Tspeed^3/48 Tspeed^2/8 - Tspeed/2 +1],[Tspeed^3/48 + Tspeed^2/8 + Tspeed/2 +1]);

OL_1DOF = F*G;
G_load  = -G/(1+F*G);

% already plotting some stuff for checking (Bode diagrams)
figure
margin(OL_1DOF)
figure
margin(OL_1DOF*Lag)
figure
margin(G_load*Lag)

OL_2DOF = (Fr+F)*G;
OL_2DOF_Delay = OL_2DOF*Lag;


%% Extra settings related to the references
ramp_rate = 5;

%% plotting some stability stuff

figure; 
nyquistplot(G,'-*'); hold on;
nyquistplot(OL_1DOF,'g-d'); hold on; 
nyquistplot(OL_2DOF,'r--'); hold on; 
nyquistplot(OL_2DOF_Delay,'k-.');  

% unit circle to be added to Nyquist
ang = 0:0.01:2*pi; 
xp  = 1*cos(ang);    
yp  = 1*sin(ang);
plot(-1+xp,0+yp,'r');
xlim([-1.5 0.1]);   ylim([-1.5 1.5])
axis equal; hold on; title('Nyquist stability')
legend('Plant','Plant+1DOF','Plant+2DOF','Plant+2DOF+Delay'); 

close(4)

%% SENSORLESS KALMAN FILTER

% noise generator seeds
seed_1 = round(rand*(2^32-1));
seed_2 = round(rand*(2^32-1));
seed_3 = round(rand*(2^32-1));
seed_4 = round(rand*(2^32-1));
seed_5 = round(rand*(2^32-1));
seed_6 = round(rand*(2^32-1));
seed_7 = round(rand*(2^32-1));
seed_8 = round(rand*(2^32-1));
seed_9 = round(rand*(2^32-1));
seed_10 = round(rand*(2^32-1));
seed_11 = round(rand*(2^32-1));

% system noise continuous standard deviation
sd_i    = 0.05; % currents standard deviation
varC_i  = sd_i^2;

% -------------------------------------------------------------------------
% digital standard deviation from Ortombina "(Automatic Tuning Procedure 
% at Standstill for Extended Kalman Filter in Sensorless Control of
% Permanent Magnet Synchronous Motors)"

% let's assume we have 50% overtorque
varD_w = 1/300*Ts^2/inertia^2*(LT*1.2)^2; % speed digital variance
varD_theta = 100/3*p^2*LT^2*Ts^4/inertia^2;

% Not really useful as a paper tbh...maybe is just my opinion...
% -------------------------------------------------------------------------

% system measurement definition (for the currents)
meas_sd = 0.1;
% init covariance matrices
Q   = eye(4);
Qd  = eye(4);
% complete the digital covariance (from symbolic)
Qd(1,1) = 0.3; %Ts*(varC_i*(R*Ts/Ld-1)^2 + Lq^2*Ts^2*1*varD_w/Ld^2);
Qd(2,2) = 0.4; %Qd(1,1); % this term has an additional one coming from the sigma_w, but it is neglected
Qd(3,3) = 5;
Qd(4,4) = 0.1;
% complete the continuous system covariance (from symbolic worked backwards)
Q(1,1) = varC_i;
Q(2,2) = varC_i+0.001;
Q(3,3) = varD_w/Ts;
Q(4,4) = (Qd(4,4) - Ts^3*Q(3,3))/Ts;

% measurement 
R_cov = meas_sd^2*eye(2);

% P0 initial covariance matrix
P0 = diag([0.1 0.1 200 10]);
%P0 = eye(4); 

%% AMSTROM APPROACH FOR ANTIWINDUP CONTROL
% the best tuning method so far found was in PID Controllers: Theory.
% Design and Tuning, by Amstrom and Hagglund. These gains are considering a
% fully digital system
% TODO: upload in the repository the .m scripts for these and not keep it
% mysterious.

kp = 5.347;
ki = (-5.135 + 5.347)/Ts;