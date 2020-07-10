function fit = RLcurve(x,level)
% function fit = RLcurve(x,level)
%
% MG Heinz
%
% Basic RL-curve model for fitting
%
% Rsp=x(1);
% Rsat=x(2);
% threshE10=x(3); % dB SPL
% threshInhib=x(4); % dB SPL
% beta=x(5); % 1-2alpha
%

% Physiological parameters
Rsp=x(1);
Rsat=x(2);
threshE10_dB=x(3); % dB SPL
threshInhib_dB=x(4); % dB SPL
beta=x(5);

% Convert physiological params into model values
Rmax=Rsat-Rsp;

thetaE=20e-6*(10.^(threshE10_dB/20))*(10/(Rmax-10))^(-1/1.77); %1/2 point
Pec=20e-6*10.^(level/20);
thetaI=20e-6*10.^(threshInhib_dB/20);
alpha=(1-beta)/2;

% MODEL
Pbm=Pec.*(1./(1+(Pec/thetaI).^2)).^alpha;
Rd=Rmax*(Pbm/thetaE).^1.77./(1+(Pbm/thetaE).^1.77);
Rtot=Rd+Rsp;

fit=Rtot;
