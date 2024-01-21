% This file is for analysing T1 for Glycerin
% Time in ms Volatge in Volts

GlycerinextractT2 = 'glycrinT21527665.isf';
 GlycerinextractT2 = isfread(GlycerinextractT2);
time_extracted = GlycerinextractT2.x;
voltage_extracted = GlycerinextractT2.y;

% This is measured data analysis
time_measured = GlycerinT2(:,1);
time_measured = table2array(time_measured);

voltage_measured =  GlycerinT2(:,2);
voltage_measured = table2array(voltage_measured);

[F,goF,Fit_output] = fit(time_measured,voltage_measured,'(a*(exp(-x/b)))+c');

err = 0.08*ones(size(voltage_measured));
errorbar(time_measured, voltage_measured, err,'vertical',"k.",'HandleVisibility','off')
hold on
plot(F,time_measured,voltage_measured)
xlabel('Time (ms)','FontSize',18);
ylabel('Potential Observed (V)','FontSize',18);
title ('Variation in Pulse Amplitude in accordance with the Nth B Pulse ','FontSize',18)
 
legend('Spin-Spin Time Delay for Glycerin', 'T2 = (30.0 +/- 3.5 )ms','FontSize',18)

% Extract weighted jacobian
J = Fit_output.Jacobian;

%Get the covariance and curvature matrix and extract the errors on F
%parameters from there.
curvature_matrix = J'*J;
covariance_matrix = inv(curvature_matrix);

% Calculate CHI_squared
min_chi2 = goF.sse;
dof = goF.dfe;

reduced_chi2 = min_chi2/dof;

figure(2)
plot(time_extracted, voltage_extracted)
xlabel('Time (s)','FontSize',18);
ylabel('Potential Observed (V)','FontSize',18);
title ('Pulse Echo signal Obtained from interaction of 15 Pulses ','FontSize',18)
 