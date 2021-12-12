function FRF = FRF_fast(x,y,fs,f_ini,f_fin)

    nn = length(x);
    nfft = 2^nextpow2(nn);

    % PODRIAMOS DEJAR LOS VALORES POR DEFECTO, PARA LA CHIRP 
    % NO HACER FALTA PERO PARA LA RANDOM SI
    window = [nn];
    noverlap = 0; %[];

    [Txy1,Ft1] = tfestimate(x,y,window,noverlap,nfft,fs);  
    
    pos1 = find(Ft1 >= f_ini); pos1 = pos1(1);
    pos2 = find(Ft1 >= f_fin); pos2 = pos2(1);

    
    % AMPLITUD Y FASE
    FRF = [Ft1(pos1:pos2),Txy1(pos1:pos2)];
    %FRF = [Ft1,Txy1];

end