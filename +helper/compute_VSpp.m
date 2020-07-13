function VSpp= compute_VSpp(spike_trains, CF_Hz, tStart, tEnd)

tStart_cell= repmat({tStart}, numel(spike_trains),1);
tEnd_cell= repmat({tEnd}, numel(spike_trains),1);
spike_trains= cellfun(@(x, y, z) make_column_vector(x(x>y & x<z)), spike_trains, tStart_cell, tEnd_cell, 'UniformOutput', false);

CF_Hz_cell= repmat({CF_Hz}, numel(spike_trains),1);

VS_trial= cell2mat(cellfun(@(x,y) get_VS_trial(x, y), spike_trains, CF_Hz_cell, 'UniformOutput', false));
Phi_trial= cell2mat(cellfun(@(x,y) get_Phi_trial(x, y), spike_trains, CF_Hz_cell, 'UniformOutput', false));
Phi_total= get_Phi_trial(cell2mat(spike_trains), CF_Hz);
VSpp= mean(VS_trial.*cos(Phi_trial-Phi_total));
end

function VS_trial= get_VS_trial(spike_train, CF_Hz)
theta_i= 2*pi*mod(spike_train, 1/CF_Hz)*CF_Hz;
VS_trial= sqrt(sum(cos(theta_i)).^2 + sum(sin(theta_i)).^2)/numel(spike_train);
end

function Phi_trial= get_Phi_trial(spike_train, CF_Hz)
theta_i= 2*pi*mod(spike_train, 1/CF_Hz)*CF_Hz;
Phi_trial= atan2(sum(sin(theta_i)), sum(cos(theta_i)));
end

function y= make_column_vector(x)
y= x(:);
end