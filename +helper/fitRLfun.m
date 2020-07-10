% function [RLVparams,exitflag]=fitRLfun(levels,rates,PlotVar)
% PlotVar=1 (or +ve) to plot the level-rate and level-rate_estimated from
% model


function [RLVparams,rates_est,exitflag]=fitRLfun(levels,rates,PlotVar,PlotFittedRLV)
if nargin==2
    PlotVar=0;
    PlotFittedRLV=0;
end

RLVparams=struct('R_SP',0,'R_Sat',50,'threshE10_dB',0,'Theta_I',0,'alpha',1/3);

x_start=[50 200  20  50 0.33];
x_lb=[  0   0   0   0   0];
x_ub=[150 400 100 100 1.0];
options=optimset('Display','none','MaxIter',1000,'MaxFunEvals',10000, 'Diagnostics','off','TolFun',1.0000e-05,'TolX',1.0000e-03);


[x_sol,~,~,exitflag,~]=lsqcurvefit(@helper.RLcurve,x_start,levels,rates,x_lb,x_ub,options);

rates_est=feval('helper.RLcurve',x_sol,levels);

if PlotVar || PlotFittedRLV
    if PlotVar
        plot(levels,rates);
        if PlotFittedRLV
            hold on;
            plot(levels,rates_est);
%             legend('Simulated Data','Model-Fit Data','Location','southeast');
        else
%         legend('Simulated Data','Location','southeast');
        end
    else
        plot(levels,rates_est);
%         legend('Model-Fit Data','Location','southeast');
    end
    
    title('Rate-Level Curve');
    xlabel('Intensity (dB SPL)');
    ylabel('Total Rate');
end

RLVparams.R_SP=x_sol(1);
RLVparams.R_Sat=x_sol(2);
RLVparams.threshE10_dB=x_sol(3);
RLVparams.Theta_I=x_sol(4);
RLVparams.alpha=(1-x_sol(5))/2;
RLVparams.levels=levels;
RLVparams.rates=rates;
RLVparams.rates_est=rates_est;

% Exitflag: 999 = not fit yet
%            -3 = COMPLEX fit
%            -2 = FIT ERROR
%            -1 = no convergence
%             0 = Max Iterations HIT
%             1 = GOOD fit
%             2 = CHECK FIT, BOUNDARY conditions HIT
%             3 = CHECK FIT, looks like params NOT ROBUST
%             8 = FIT ==> RL model is not appropriate

if(~isreal(x_sol))
    exitflag=-3;
    disp('COMPLEX SOLUTION');
elseif(exitflag==-2)
    disp('ERROR IN FITTING');
elseif(exitflag==-1)
    disp('DID NOT CONVERGE');
elseif(exitflag==0)
    disp('MAX ITERATIONS MADE');
end

for i=1:5
    if((round(100*x_sol(i))==round(100*x_lb(i)))||(round(100* x_sol(i))==round(100*x_ub(i))))
        exitflag=2;
%         fprintf('***Boundary condition restricted Solution (for param %d)***',i);
    end
end