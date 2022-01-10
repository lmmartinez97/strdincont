clearvars; 
close all
clc

set(0,'defaultAxesFontName',      'Times New Roman')
set(0,'defaultTextFontName',      'Times New Roman')
set(0,'defaultUicontrolFontName', 'Times New Roman')
set(0,'defaultUitableFontName',   'Times New Roman')
set(0,'defaultAxesFontName',      'Times New Roman')
set(0,'defaultTextFontName',      'Times New Roman')
set(0,'defaultUipanelFontName',   'Times New Roman')

tam = 12;  
set(0, 'DefaultAxesFontSize',tam)
set(0, 'DefaultTextFontSize',tam)

% Structure data
f1 = 1.4;                     
damp1 = 2/100;            
m1 = 2000; 
k1 = m1*((f1*2*pi)^2);
c1 = m1*2*damp1*(f1*2*pi);

u_active = 200;

% CHIRP EXCITATION
fmin = 0.1;
fmax = 2.5;   
tend = 500;

fs = 1500;
ts = 1/fs; 
tsam = ts;

% TMD data
mu = 0.03;  % mass ratio
mt = mu*m1; % TMD mass
sat = 50; %saturates at about 50 N
dampt = sqrt(3*mu/(8*(1+mu)))*sqrt(1+mu*27/32);
ft = f1*sqrt(1/(1+mu));               

kt = ((2*pi*ft)^2)*mt;                          
ct = 2*mt*(2*pi*ft)*dampt; 
I1 = 0.01;

% The Structure change (sólo si queremos)  
f1d = f1*1;               
damp1d = damp1;      
m1d = m1; 
k1d = m1d*((f1d*2*pi)^2);
c1d = m1d*2*damp1d*(f1d*2*pi);

% Transfer function Structure without TMD
s = tf('s');
Gsd = ((1/m1d)*s*s)/(s*s+s*(c1d/m1d) + (k1d/m1d));
[N_Gs,D_Gs] = tfdata(Gsd,'v');

% To semi-active control law
Kmax = 1*ct;
Kmin = ct/2;

% EXECUTE MODEL
sim control_STMD_Ideal_Entrega_Final.slx

% Plot results
% Time histories
mm(1) = max(abs(acc));
mm(2) = max(abs(acc_TMD));
mm(3) = max(abs(acc_STMD));
mm(4) = max(abs(acc_STMD1));

figure
    X = categorical({'uncontrolled', 'w TMD', 'w STMD', 'w STMD1'});
    X = reordercats(X,{'uncontrolled', 'w TMD', 'w STMD', 'w STMD1'});
    bar(X, mm)
    xlabel('Controller type')
    ylabel('Max acceleration [m/s^2]')
    title('Acceleration comparison')

figure
subplot(2,1,1)
    plot(t,acc,'Color',[0.8 0.8 0.8])
    hold on
    plot(t,acc_TMD,'r')
    plot(t,acc_STMD,'b')
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    legend({'Uncontrolled','w TMD','w STMD'},'Location','Southeast')

subplot(2,1,2)
    plot(t,acc,'Color',[0.8 0.8 0.8])
    hold on
    plot(t,acc_TMD,'r')
    plot(t,acc_STMD1,'m')
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    legend({'Uncontrolled','w TMD', 'w STMD1'},'Location','Southeast')

% Compute experimental FRFs
FRF1 = FRF_fast(input,acc,fs,fmin,fmax);
FRF2 = FRF_fast(input,acc_TMD,fs,fmin,fmax);
FRF3 = FRF_fast(input,acc_STMD,fs,fmin,fmax);
FRF4 = FRF_fast(input,acc_STMD1,fs,fmin,fmax);

figure
subplot(2,1,1)
    plot(FRF1(:,1),20*log10(abs(FRF1(:,2))),'Color',[0.8 0.8 0.8],'Linewidth',1)
    hold on
    plot(FRF2(:,1),20*log10(abs(FRF2(:,2))),'r','Linewidth',1)
    plot(FRF3(:,1),20*log10(abs(FRF3(:,2))),'b','Linewidth',1)
    plot(FRF4(:,1),20*log10(abs(FRF4(:,2))),'m','Linewidth',1)

    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]') 
    xlim([f1d-0.2 f1d+0.2])
    legend({'Uncontrolled','w TMD','w STMD', 'w STMD1'},'Location','Southeast')

subplot(2,1,2)
    plot(FRF1(:,1),unwrap((angle(FRF1(:,2))),10).*180/pi,'Color',[0.8 0.8 0.8],'Linewidth',1)
    hold on
    plot(FRF2(:,1),unwrap((angle(FRF2(:,2))),10).*180/pi,'r','Linewidth',1)
    plot(FRF3(:,1),unwrap((angle(FRF3(:,2))),10).*180/pi,'b','Linewidth',1)
    plot(FRF4(:,1),unwrap((angle(FRF4(:,2))),10).*180/pi,'m','Linewidth',1)
    line([f1d-0.2 f1d+0.2],[90 90],'Linestyle','-.','Color',[0.5 0.5 0.5])
    xlabel('Frequency [Hz]')
    ylabel('Phase [º]') 
    xlim([f1d-0.2 f1d+0.2])
    legend({'Uncontrolled','w TMD','w STMD', 'w STMD1'},'Location','Southeast')


disp("Press any key to continue with structure modification")
pause;

% The Structure change (sólo si queremos)  
f1d = f1*0.9;               
damp1d = damp1;      
m1d = m1; 
k1d = m1d*((f1d*2*pi)^2);
c1d = m1d*2*damp1d*(f1d*2*pi);

% Transfer function Structure without TMD
s = tf('s');
Gsd = ((1/m1d)*s*s)/(s*s+s*(c1d/m1d) + (k1d/m1d));
[N_Gs,D_Gs] = tfdata(Gsd,'v');

% EXECUTE MODEL
sim control_STMD_Ideal_Entrega_Final.slx

% Plot results
% Time histories
mm(1) = max(abs(acc));
mm(2) = max(abs(acc_TMD));
mm(3) = max(abs(acc_STMD));
mm(4) = max(abs(acc_STMD1));

figure
    X = categorical({'uncontrolled', 'w TMD', 'w STMD', 'w STMD1'});
    X = reordercats(X,{'uncontrolled', 'w TMD', 'w STMD', 'w STMD1'});
    bar(X, mm)
    xlabel('Controller type')
    ylabel('Max acceleration [m/s^2]')
    title('Acceleration comparison')

figure
subplot(2,1,1)
    plot(t,acc,'Color',[0.8 0.8 0.8])
    hold on
    plot(t,acc_TMD,'r')
    plot(t,acc_STMD,'b')
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    legend({'Uncontrolled','w TMD','w STMD'},'Location','Southeast')

subplot(2,1,2)
    plot(t,acc,'Color',[0.8 0.8 0.8])
    hold on
    plot(t,acc_TMD,'r')
    plot(t,acc_STMD1,'m')
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    legend({'Uncontrolled','w TMD', 'w STMD1'},'Location','Southeast')


% Compute experimental FRFs
FRF1 = FRF_fast(input,acc,fs,fmin,fmax);
FRF2 = FRF_fast(input,acc_TMD,fs,fmin,fmax);
FRF3 = FRF_fast(input,acc_STMD,fs,fmin,fmax);
FRF4 = FRF_fast(input,acc_STMD1,fs,fmin,fmax);

figure
subplot(2,1,1)
    plot(FRF1(:,1),20*log10(abs(FRF1(:,2))),'Color',[0.8 0.8 0.8],'Linewidth',1)
    hold on
    plot(FRF2(:,1),20*log10(abs(FRF2(:,2))),'r','Linewidth',1)
    plot(FRF3(:,1),20*log10(abs(FRF3(:,2))),'b','Linewidth',1)
    plot(FRF4(:,1),20*log10(abs(FRF4(:,2))),'m','Linewidth',1)

    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]') 
    xlim([f1d-0.2 f1d+0.2])
    legend({'Uncontrolled','w TMD','w STMD', 'w STMD1'},'Location','Southeast')

subplot(2,1,2)
    plot(FRF1(:,1),unwrap((angle(FRF1(:,2))),10).*180/pi,'Color',[0.8 0.8 0.8],'Linewidth',1)
    hold on
    plot(FRF2(:,1),unwrap((angle(FRF2(:,2))),10).*180/pi,'r','Linewidth',1)
    plot(FRF3(:,1),unwrap((angle(FRF3(:,2))),10).*180/pi,'b','Linewidth',1)
    plot(FRF4(:,1),unwrap((angle(FRF4(:,2))),10).*180/pi,'m','Linewidth',1)
    line([f1d-0.2 f1d+0.2],[90 90],'Linestyle','-.','Color',[0.5 0.5 0.5])
    xlabel('Frequency [Hz]')
    ylabel('Phase [º]') 
    xlim([f1d-0.2 f1d+0.2])
    legend({'Uncontrolled','w TMD','w STMD', 'w STMD1'},'Location','Southeast')


disp('Press enter to continue with sinusoidal simulations')
pause;

for i=1:5    freq = (2-0.1)/5 * i;
    fmin = freq;
    fmax = freq;
    tend = 15;
    % The Structure change (sólo si queremos)  
    f1d = f1*1;               
    damp1d = damp1;      
    m1d = m1; 
    k1d = m1d*((f1d*2*pi)^2);
    c1d = m1d*2*damp1d*(f1d*2*pi);
    
    % Transfer function Structure without TMD
    s = tf('s');
    Gsd = ((1/m1d)*s*s)/(s*s+s*(c1d/m1d) + (k1d/m1d));
    [N_Gs,D_Gs] = tfdata(Gsd,'v');

    sim control_STMD_Ideal_Entrega_Final.slx
    legendCell{i} = strcat('f = ', num2str(freq));
    figure(15)
    subplot(2,1,1)
        hold on
        plot(t,acc_STMD)
        xlabel('Time [s]')
        ylabel('Acceleration [m/s^2]')
    
    subplot(2,1,2)
        hold on
        plot(t,acc_STMD1)
        xlabel('Time [s]')
        ylabel('Acceleration [m/s^2]')
end
figure(15)
    subplot(2,1,1)
    legend(legendCell)
    subplot(2,1,2)
    legend(legendCell)
