% Plot results
% Time histories
max(abs(out.acc));
%max(abs(out.acc_TMD));
max(abs(out.acc_STMD));
figure
    plot(t,out.acc,'-r')
    hold on
    %plot(t,acc_TMD,'r')
    plot(t,out.acc_STMD,'.-b')
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    legend({'Uncontrolled','w STMD'},'Location','Southeast')

% Compute experimental FRFs
FRF1 = FRF_fast(out.input,out.acc,fs,fmin,fmax);
%FRF2 = FRF_fast(input,acc_TMD,fs,fmin,fmax);
FRF3 = FRF_fast(out.input,out.acc_STMD,fs,fmin,fmax);

figure
    plot(FRF1(:,1),20*log10(abs(FRF1(:,2))),'Color',[0.8 0.8 0.8],'Linewidth',1)
    hold on
    %plot(FRF2(:,1),20*log10(abs(FRF2(:,2))),'r','Linewidth',1)
    plot(FRF3(:,1),20*log10(abs(FRF3(:,2))),'b','Linewidth',1)
    xlabel('Frequency [Hz]')
    ylabel('Magnitude [dB]') 
    xlim([f1d-0.5 f1d+0.5])
    legend({'Uncontrolled','w STMD'},'Location','Southeast')

figure
    plot(FRF1(:,1),unwrap((angle(FRF1(:,2))),10).*180/pi,'Color',[0.8 0.8 0.8],'Linewidth',1)
    hold on
    %plot(FRF2(:,1),unwrap((angle(FRF2(:,2))),10).*180/pi,'r','Linewidth',1)
    plot(FRF3(:,1),unwrap((angle(FRF3(:,2))),10).*180/pi,'b','Linewidth',1)
    line([0.6 1.5],[90 90],'Linestyle','-.','Color',[0.5 0.5 0.5])
    xlabel('Frequency [Hz]')
    ylabel('Phase [º]') 
    xlim([f1d-0.2 f1d+0.2])
    legend({'Uncontrolled','w STMD'},'Location','Southwest')