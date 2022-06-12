 for M1=10:10:40
    NFFT = 1024;
    [PBT, fgrid] = btmethod(x, M1, NFFT);
    figure(5)
    plot(fgrid,PBT)
    hold on
 end

legend('BT-M1 = 10','BT-M1 = 20','BT-M1 = 30','BT-M1 = 40','Bartlett')