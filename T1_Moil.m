% This file is for analysing T1 for Mineral Oil

% Time is in ms, Voltage is in volts

time = MoilT1data(:,1);
time = table2array(time);
err_time = (2E-08)*ones(size(time));

voltage =  MoilT1data(:,2);
voltage = table2array(voltage);
err = 0.08*ones(size(voltage));

%F = fit(time,voltage,'(a*(1-exp(-x/b)))+c','Start', [100, 2, 0],'Weight',1/0.02^2+zeros(length(time),1));
[F,goF,Fit_output] = fit(time,voltage,'(a*(1-exp(-x/b)))+c','Start', [100, 2, 0],'Weight',err.^(-2));

plot(F,time,voltage)
xlabel('Time (ms)','FontSize',18);
ylabel('Potential Observed (V)','FontSize',18);
title ('Variation in Pulse Amplitude with increasing time Delay between Pulses','FontSize',18)

legend('Spin-Lattice Time Delay for Mineral Oil', 'T1 = (22.45 +/- 0.16)ms','FontSize',18)
hold on
errorbar(time, voltage, err,"k.",'HandleVisibility','off')
hold on
errorbar(time, voltage, err_time,'horizontal', "k.",'HandleVisibility','off')
hold off

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

