function [lags,Twindow] = DFC_timeseg(twi,tws,over,fs,time)
    cont    = 1;
    lags    = cell(1,numel((twi:0.1:tws)));
    Twindow = cell(1,numel((twi:0.1:tws)));
    for tw = (twi:0.1:tws)
        Twindow{1,cont} = round(tw*fs);
        overlap         = round(over*Twindow{1,cont});
        lags{1,cont}    = 0:overlap:(time*fs)-Twindow{1,cont};
        cont            = cont + 1;
    end
end