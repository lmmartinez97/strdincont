clearvars; 
close all
clc

% Structure data
f1 = 1;                     
damp1 = (1)/100;            
m1 = 1; 
k1 = m1*((f1*2*pi)^2);
c1 = m1*2*damp1*(f1*2*pi);

% CHIRP EXCITATION
fmin = 0.001;
fmax = 2;   
tend = 500;

fs = 1500;
ts = 1/fs; 

% TMD data
mu = 0.02;  % mass ratio
mt = mu*m1; % TMD mass

dampt = sqrt(3*mu/(8*(1+mu)))*sqrt(1+mu*27/32);
ft = f1*sqrt(1/(1+mu));               

kt = ((2*pi*ft)^2)*mt;                          
ct = 2*mt*(2*pi*ft)*dampt; 
sat = 50;

% The Structure change (sólo si queremos)  
f1d = f1*1.2;               
damp1d = damp1;      
m1d = m1; 
k1d = m1d*((f1d*2*pi)^2);
c1d = m1d*2*damp1d*(f1d*2*pi);

% Transfer function Structure without TMD
s = tf('s');
Gsd = ((1/m1d)*s*s)/(s*s+s*(c1d/m1d) + (k1d/m1d));
[N_Gs,D_Gs] = tfdata(Gsd,'v');

% To semi-active control law
Kmax = 50*ct;
Kmin = ct/2;
