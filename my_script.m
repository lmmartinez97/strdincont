clearvars; 
close all
clc

% Structure data
fs = 1;                     
damps = (1)/100;            
ms = 1; 
ks = ms*((fs*2*pi)^2);
cs = ms*2*damps*(fs*2*pi);

% CHIRP EXCITATION
fmin = 0.001;
fmax = 2;   
tend = 500;

fs = 1500;
ts = 1/fs; 

% TMD data
mu = 0.02;  % mass ratio
md = mu*ms; % TMD mass

dampt = sqrt(3*mu/(8*(1+mu)))*sqrt(1+mu*27/32);
fd = fs*sqrt(1/(1+mu));               

kd = ((2*pi*fd)^2)*md;                          
cd = 2*md*(2*pi*fd)*dampt; 
sat = 50;

% The Structure change (sólo si queremos)  
fsd = fs*1.2;               
dampsd = damps;      
msd = ms; 
ksd = msd*((fsd*2*pi)^2);
csd = msd*2*dampsd*(fsd*2*pi);

% Transfer function Structure without TMD
s = tf('s');
Gsd = ((1/msd)*s*s)/(s*s+s*(csd/msd) + (ksd/msd));
[N_Gs,D_Gs] = tfdata(Gsd,'v');

% To semi-active control law
Kmax = 50*cs;
Kmin = cs/2;


A = [      0          1       0    0      0    0;
           0          0       1    0      0    0;
     -(ks+kd)/ms  -(cs-cd)/ms 0  kd/ms  cd/ms  0;
           0          0       0    0      1    0;
           0          0       0    0      0    1;
         kd/md      cd/ms     0  kd/md  cd/md  0
     ];

B = [0 0 1/ms 0 0 0]';

C = eye(6);

D = zeros(6,1);
