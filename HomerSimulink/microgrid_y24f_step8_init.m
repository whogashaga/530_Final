%% Microgrid Simulation
% Created: ECE 530 class of fall 2024

clc
close all   % close figure windows
clear
format compact



%% Simulation settings

simu.endTime = 60*60*24;
simu.maxStepSize = 1e-1;



%% Data

load windspeedtimeseries
load loadtimeseries
load illuminationcurrenttimeseries

meanWindSpeed = 8;
maxWindSpeed = 15;

meanLoadPower = 1.98e3; %Replace 80e3 with 1.98e3 in your code to ensure the units are consistent with "kWh/day = 47.51.
maxLoadPower = 2.36e3; %Replace 100e3 with 2.36e3 in your code to match the peak load of 56.60 kWh/day.

maxIlluminationCurrent = 10; % 10 is fully sunny day

% Scaling
windspeedtimeseries.Data = fn_scaling(...
    meanWindSpeed*windspeedtimeseries.Data,meanWindSpeed,maxWindSpeed);
loadtimeseries.Data = fn_scaling(...
    loadtimeseries.Data,meanLoadPower,maxLoadPower);
illuminationcurrenttimeseries.Data = ...
    illuminationcurrenttimeseries.Data*maxIlluminationCurrent/max(illuminationcurrenttimeseries.Data);



%% Wind Turbine Parameters

wt.N_turbines = 0; % We don't use it

wt.rho = 1.2; % Density of air kg/m^3

% XANT 100 kW
wt.Pgen_rated = 100e3;
wt.Cfric = 0.01*100e3/6^2;    % Nm/(rad/s)    B*6^2 = 0.01*100e3
wt.bladeLength = 11;
wt.bladeWeight = 1000;    % From http://windpower.sandia.gov/other/041605C.pdf
wt.J = 3*1/3*wt.bladeWeight*wt.bladeLength^2;  % Missing generator, gearbox, shaft, etc...

wt.A = pi*wt.bladeLength^2;

wt.w_0 = 0;

% Cp curve modeling
% lambdaai = 1/( (1/(lambda+0.08*beta)) - 0.035/(beta^3+1) )
% cp = c1*(c2/lambdaai-c3*beta-c4)*exp(-c5/lambdaai) + c6*lambda
plot_cp = 0;
if plot_cp
    wt.c = [
        0.5176
        116
        0.4
        5
        21
        0.0068
        ];
    figure
    lambda = [0:0.1:13];
    for beta = [0:5:30];
        lambdaai = 1./(1./(lambda+0.08*beta) - 0.035./(beta.^3+1));
        cp = wt.c(1)*(wt.c(2)./lambdaai - wt.c(3)*beta - wt.c(4)) .* exp(-wt.c(5)./lambdaai) + wt.c(6)*lambda;
        hold on
        plot(lambda,cp)
        hold off
    end
    axis([0 13 0 0.5])
    xlabel('lambda (tip speed ratio)')
    ylabel('Cp')
end

% From Cp(lambda) plot
wt.lambda_opt = 8.1;
wt.Cp_opt = 0.48;

% Region 2 and 3 boundary -> rated rotational speed and wind speed
wt.u_rated = (wt.Pgen_rated/(wt.Cp_opt*0.5*1.2*wt.A))^(1/3);   % P_rated = Cp_max*0.5*A*bladeLength*u_rated^3
wt.w_rated = wt.lambda_opt*wt.u_rated/wt.bladeLength;
wt.Tgen_rated = wt.Pgen_rated/wt.w_rated;

wt.pitchctrl.kp = 1;
wt.pitchctrl.ki = 0.1;
wt.pitchctrl.upperLimit = 1;
wt.pitchctrl.lowerLimit = 0;



%% Energy Storage

es.eta_pe = 0.95;
es.eta_sm = sqrt(0.9);
es.SOC_0 = 1; % Updated SOC_0 to 1.0 to reflect a 100% initial state of charge.

es.P_pe_rated = 120e3;
es.P_pe_slew_upper = es.P_pe_rated / 1;
es.P_pe_slew_lower = -es.P_pe_rated / 1;

es.E_rated_kWh = 1; % Adjusted E_rated_kWh to 1.0 for the battery's rated energy capacity (1 kWh).

es.E_rated = es.E_rated_kWh*1000*3600;

es.P_rated = es.P_pe_rated;



%% Microgrid

mg.H = 5;   
mg.D = 1;
mg.P_base = 2.36e3; %Replace 100e3 with 2.36e3 in your code to match the peak load of 56.60 kWh/day.
mg.wpu_0 = 1;
mg.X = 0.05;

es.Kgpri = 20;

ht.Kgpri = 20;
ht.Kgsec = 20/60;



%% Solar

% Parameters based on SolarWorld 300 module
% 300 W panel, 60 cells, 5 W per cell
% Open circuit voltage of 40 V (0.66 V per cell)
% Short circuit current of 9.83 A
% Maximum power point at 32.6 V (0.54 V per cell) and 9.31 A

pv.P_rated = 10e3; % Updated to reflect the rated capacity of 10 kW

pv.Is = 1e-10; % Produces about 0.66 V at 9.8 A of current
pv.Rs = 0.005; % Adjusted this to get the approximate max power point
pv.Rp = 2500;  % Based on L. Ma et al. "The Measurement of Series and Shunt Resistances of the Silicon Solar Cell Based on LabVIEW"
pv.VT = 0.026;

pv.vd_0 = 0.7;

pv.MPPT_sampleTime = 1;  % Not optimized

pv.N_photocells = 60;
pv.N_panels = pv.P_rated/300;



%% Hydro

ht.P_rated = 100e3;

ht.tauw = 5;
ht.q_0 = 0.1;
ht.beta = 10;

ht.powercontroller.kp = 1/4;         
ht.powercontroller.ki = ht.powercontroller.kp/2;  
ht.powercontroller.int_0 = 0.1;
ht.powercontroller.kt = 1;
ht.powercontroller.lowerLimit = 0;    
ht.powercontroller.upperLimit = 1;  

ht.servo.tau = 0.1;
ht.servo.speedUpperLimit = 0.1;
ht.servo.speedLowerLimit = -0.1;
ht.servo.posUpperLimit = 1;
ht.servo.posLowerLimit = 0.05;
ht.servo.initPos = 0.1;



%% Water Pump (Deferrable Load) 
% We don't need to use it
wp.P_rated = meanLoadPower;
wp.height = 40;

wp.P_upper = wp.P_rated;
wp.P_lower = 0;
wp.P_upperRate = wp.P_rated/10;
wp.P_lowerRate = -wp.P_upperRate;



%% Closing

disp("Good job!  The init file ran succesfully, hopefully the simulation does too.")