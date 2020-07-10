function spl_out= calc_dbspl(vecin)

pRef= 20e-6;
vecin= vecin-mean(vecin);
spl_out= 20*log10(rms(vecin)/pRef);