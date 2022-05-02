% Implementation of Sliding Mode Observer sensorless according to the papers
% given in the following paper: 
% "Extended Range PMSM Sensorless Speed Drive Based on Stochastic Filtering"

% by Nicola Chiodetto @02/05/2022

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

% Controller modeling -----------------------------------------------------
% SPEED CONTROLLER --------------------------------------------------------
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

%% AMSTROM APPROACH FOR ANTIWINDUP CONTROL
% the best tuning method so far found was in PID Controllers: Theory.
% Design and Tuning, by Amstrom and Hagglund. These gains are considering a
% fully digital system
% TODO: upload in the repository the .m scripts for these and not keep it
% mysterious.

kp = 5.347;
ki = (-5.135 + 5.347)/Ts;

% Sliding mode observer parameters (Refer to the multitude of SMO papers)--
% estimator proportionality factors
f       = -R/Ld;
g       = 1/Ld;

% PLL ---------------------------------------------------------------------
% Gain definition ---------------------------------------------------------
zeta    = 1;
alpha   = 0.2; 
w_n     = 250*2*pi*alpha; 
Ki_est  = w_n^2;
Kp_est  = 2*zeta*w_n;

% speed filter
speedFilter_freq = 75;

% comparing the tf for open loop of speed, pll and filter of first order
speed_tf    = feedback(OL_1DOF,1);
PLL_tf      = feedback(tf([Kp_est Ki_est],[1 0 0]),1);
filter_tf   = tf([2*pi*speedFilter_freq],[1 2*pi*speedFilter_freq]);

figure
bode(speed_tf);     hold on
bode(PLL_tf);       hold on
bode(filter_tf);    hold off
legend('speed loop','pll','filter')