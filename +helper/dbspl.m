% Important: to calculate dbspl of a signal, use calc_dbspl
% This function converts linear amplitude to dbSPL
function vec_out= dbspl(vec_in)


pRef= 20e-6;
vec_out= 20*log10(vec_in/pRef);