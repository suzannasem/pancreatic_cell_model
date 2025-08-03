
% This file runs the complete alpha-, beta-cell, and delta-cell model in a perfusion
% setting with different step changes in glucose in flow rates.

%Experimental conditions
numIslets = 5; %Number of islets
g_in_0 = 1; %mM - Initial glucose in flow concentration
betaOnOff = true; %true = on, false = off
alphaOnOff = true; %true = on, false = off
deltaOnOff = true; 
I_t_in_c = 0; %mg/dL - Insulin in flow step concentration
G_t_in_c = 0; %mg/dL - Glucagon in flow step concentration
S_t_in_c = 0; % mg/dL - somatostatin in flow step conc

g_up_vals =  (1:0.1:180)'; %mM - range of glucose values to examine


%Parameters
%--------------------------------------------------------------------------
%Basal levels for humans
gba_ = 90.08; %mg/dL - basal glucose levels (Alcazar & Buchwald, 2019)
Gba_ = 2.5e-6; %mg/dL - basal glucagon levels (van Vliet et. al., 2020 and http://www.ncbi.nlm.nih.gov/books/NBK279127)
Iba_ = 2.2e-5; %mg/dL - basal insulin levels (Van Vliet et. al., 2020)

% NEW
Sba_ = 1.3e-6; % mg/dL - basal somatostatin levels (Conlon 1983) https://diabetesjournals.org/diabetes/article/32/8/723/5736/Circulating-Somatostatin-Concentrations-in-Healthy
% goal: plot glucagon varying somatostatin basal levels
Sba2_ = 2*1.3e-6; % 2x
Sba3_ = 3*1.3e-6; %3x
Sba4_ = 0.5*1.3e-6; % 0.5x
%--------------------------------------------------------------------------
%Kinetic secretion parameters
k_gB_ = 0.554; %1/min - beta cell glucose transduction
k_G_  = 0.554; %1/min - beta cell glucagon transduction
k_gA_ = 0.554/25; %1/min - alpha cell glucose transduction
k_I_  = 0.554*5; %1/min - alpha cell insulin transduction
k_S_  = 0.554; % 1/min - delta cell somatostatin transduction

%Rate constant for transfer from insulin pool 1 to insulin pool 2
m_I1_ = 0.336; %1/min - max value
h_I1_ = 3.75; %[] - half-maximum value
n_I1_ = 9.97; %[] - Hill exponent

%Rate constant for transfer from insulin pool 2 to outside cell
m_I2_ = 0.360; %1/min - max value
h_I2_ = 0.968; %[] - half-maximum value
n_I2_ = 6.68; %[] - Hill exponent

%Rate constant for transfer from glucagon pool 1 to glucagon pool 2
m_G1_ = 0.336; %1/min - max value
h_G1_ = 3.75; %[] - half-maximum value
n_G1_ = 9.97; %[] - Hill exponent

%Rate constant for transfer from glucagon pool 2 to outside cell
m_G2_ = 0.360; %1/min - max value
h_G2_ = 0.968; %[] - half-maximum value
n_G2_ = 6.68; %[] - Hill exponent

% NEW

%Rate constant for transfer from somatostatin pool 1 to somatostatin pool 2
%- i.e. k5
m_S1_ = 0.336; %1/min - max value
h_S1_ = 3.75; %[] - half-maximum value
n_S1_ = 9.97; %[] - Hill exponent

%Rate constant for transfer from somatostatin pool 2 to outside cell - i.e.
%k6
m_S2_ = 0.360; %1/min - max value
h_S2_ = 0.968; %[] - half-maximum value
n_S2_ = 6.68; %[] - Hill exponent

%--------------------------------------------------------------------------
%Steady-state model parameters

%Interaction parameters from mice

%Net beta-cell signal
m_GB_ = 1.11; %[] - maximum contribution of glucagon to beta cell signal
h_GB_ = 502; %[] - half-maximum value for contribution of glucagon
n_GB_ = 0.628; %[] - Hill exponent for contribution of insulin
h_gB_ = 1.07; %[] - half-maximum value for on/off effect of glucose
n_gB_ = 0.350; %[] - Hill exponent for on/off effect of glucose
n_SB_ = 0.35; % hill exponent for effect of somatostatin
h_SB_ = 1.07; % half maximal contribution of somatostatin


%Net alpha-cell signal
h_IA_ = 10.0; %[] - half-maximum for contribution of insulin to alpha cell signal
n_IA_ = 1.17; %[] - Hill exponent for contribution of insulin
m_g_ = 0.600; %[] Maximum fraction of glucose signal that insulin can remove
n_SA_ = 0.35 ; % hill coefficient of effect of somatostatin
h_SA_ = 1.07; % half maximal somatostatin contribution

%--------------------------------------------------------------------------

%Mass balance related parameters
Q_ =  0.05/10^3*10; %dL/min - perfusion rate
V_P_ = 1/10^3*10; %dL - volume

%--------------------------------------------------------------------------
%Steady-state secretion parameters

%Steady-state insulin secretion
m_I_ = 1.03e-7./15.*numIslets.*betaOnOff; %mg/min - maximum steady-state insulin secretion rate
h_I_ = 3.97; %[] - half-maximum for insulin secretion rate
n_I_ = 4.84; %[] - Hill exponent for insulin secretion rate

%Steady-state glucagon secretion
m_G_ = 2.24e-9./15.*numIslets.*alphaOnOff; %mg/min - maximum steady-state glucagon secretion rate
h_G_ = 1.06; %[] - half-maximum for glucagon secretion rate
n_G_ = 3.5; %[] - Hill exponent for glucagon secretion rate

% Steady-state somatostatin secretion
m_S_ = 1.03e-7./15.*numIslets.*deltaOnOff; %mg/min - maximum steady-state somatostatin secretion rate
h_S_ = 3.97; %[] - half-maximum for insulin secretion rate
n_S_ = 4.84; %[] - Hill exponent for insulin secretion rate

%Background signals
X_B0_ = 2.60;
X_A0_ = 4.40;

%--------------------------------------------------------------------------


params_ = [gba_, Gba_, Iba_, Sba_...
          k_gB_, k_G_, k_gA_, k_I_, k_S_...
          m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
          h_IA_, n_IA_, X_A0_, m_g_, ...
          m_I_, h_I_, n_I_, ...
          m_G_, h_G_, n_G_, ...
          m_S_, h_S_,n_S_, ...
          m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
          m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
          m_S1_, h_S1_, n_S1_, m_S2_, h_S2_, n_S2_, ...
          n_SB_,h_SB_,n_SA_,h_SA_,...
          Q_,V_P_];

% 2x
params2_ = [gba_, Gba_, Iba_, Sba2_...
          k_gB_, k_G_, k_gA_, k_I_, k_S_...
          m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
          h_IA_, n_IA_, X_A0_, m_g_, ...
          m_I_, h_I_, n_I_, ...
          m_G_, h_G_, n_G_, ...
          m_S_, h_S_,n_S_, ...
          m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
          m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
          m_S1_, h_S1_, n_S1_, m_S2_, h_S2_, n_S2_, ...
          n_SB_,h_SB_,n_SA_,h_SA_,...
          Q_,V_P_];

% 3x
params3_ = [gba_, Gba_, Iba_, Sba3_...
          k_gB_, k_G_, k_gA_, k_I_, k_S_...
          m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
          h_IA_, n_IA_, X_A0_, m_g_, ...
          m_I_, h_I_, n_I_, ...
          m_G_, h_G_, n_G_, ...
          m_S_, h_S_,n_S_, ...
          m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
          m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
          m_S1_, h_S1_, n_S1_, m_S2_, h_S2_, n_S2_, ...
          n_SB_,h_SB_,n_SA_,h_SA_,...
          Q_,V_P_];

% 0.5x
params4_ = [gba_, Gba_, Iba_, Sba4_...
          k_gB_, k_G_, k_gA_, k_I_, k_S_...
          m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
          h_IA_, n_IA_, X_A0_, m_g_, ...
          m_I_, h_I_, n_I_, ...
          m_G_, h_G_, n_G_, ...
          m_S_, h_S_,n_S_, ...
          m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
          m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
          m_S1_, h_S1_, n_S1_, m_S2_, h_S2_, n_S2_, ...
          n_SB_,h_SB_,n_SA_,h_SA_,...
          Q_,V_P_];
      
%Step functions are used, because this better represents what is done
%experimentally
G_t_in = @(t) (t > 0)*G_t_in_c; %Glucagon in flow concentration step function
I_t_in = @(t) (t > 0)*I_t_in_c; %Insulin in flow concentration step function
S_t_in = @(t) (t > 0)*S_t_in_c; % somatostatin in flow concentration step function
t = 0:0.1:60; %min - time range to examine

n = length(g_up_vals);

R_I = zeros(length(g_up_vals),1);
R_G = zeros(length(g_up_vals),1);
R_S = zeros(length(g_up_vals),1);

R_I2 = zeros(length(g_up_vals),1);
R_G2 = zeros(length(g_up_vals),1);
R_S2 = zeros(length(g_up_vals),1);

R_I3 = zeros(length(g_up_vals),1);
R_G3 = zeros(length(g_up_vals),1);
R_S3 = zeros(length(g_up_vals),1);

R_I4 = zeros(length(g_up_vals),1);
R_G4 = zeros(length(g_up_vals),1);
R_S4 = zeros(length(g_up_vals),1);

for i = 1:n

    g_t_in = @(t) g_in_0.*18.016 + (t > 0).*(g_up_vals(i) - g_in_0).*18.016;
    [~,results] = simulate_alphaBetaModel_perfusion_delta(params_,g_t_in,G_t_in,I_t_in,S_t_in,t);
    
    R_I(i) = (results(end,1) - results(1,1)).*V_P_ - Q_.*trapz(t,I_t_in(t)) + Q_.*trapz(t,results(:,1));
        %Because I_t_in is a function of time, need to calculate insulin
        %secretion by integrating the mass balance
    R_G(i) = (results(end,2) - results(1,2)).*V_P_ - Q_.*trapz(t,G_t_in(t)) + Q_.*trapz(t,results(:,2));
        %Because G_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance

    R_S(i) = (results(end,3) - results(1,3)).*V_P_ - Q_.*trapz(t,S_t_in(t)) + Q_.*trapz(t,results(:,3));
        %Because S_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance 
end

for i = 1:n

    g_t_in = @(t) g_in_0.*18.016 + (t > 0).*(g_up_vals(i) - g_in_0).*18.016;
    [~,results2] = simulate_alphaBetaModel_perfusion_delta(params2_,g_t_in,G_t_in,I_t_in,S_t_in,t);
    
    R_I2(i) = (results2(end,1) - results2(1,1)).*V_P_ - Q_.*trapz(t,I_t_in(t)) + Q_.*trapz(t,results2(:,1));
        %Because I_t_in is a function of time, need to calculate insulin
        %secretion by integrating the mass balance
    R_G2(i) = (results2(end,2) - results2(1,2)).*V_P_ - Q_.*trapz(t,G_t_in(t)) + Q_.*trapz(t,results2(:,2));
        %Because G_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance

    R_S2(i) = (results2(end,3) - results2(1,3)).*V_P_ - Q_.*trapz(t,S_t_in(t)) + Q_.*trapz(t,results2(:,3));
        %Because S_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance 
end

for i = 1:n

    g_t_in = @(t) g_in_0.*18.016 + (t > 0).*(g_up_vals(i) - g_in_0).*18.016;
    [~,results3] = simulate_alphaBetaModel_perfusion_delta(params3_,g_t_in,G_t_in,I_t_in,S_t_in,t);
    
    R_I3(i) = (results3(end,1) - results3(1,1)).*V_P_ - Q_.*trapz(t,I_t_in(t)) + Q_.*trapz(t,results3(:,1));
        %Because I_t_in is a function of time, need to calculate insulin
        %secretion by integrating the mass balance
    R_G3(i) = (results3(end,2) - results3(1,2)).*V_P_ - Q_.*trapz(t,G_t_in(t)) + Q_.*trapz(t,results3(:,2));
        %Because G_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance

    R_S3(i) = (results3(end,3) - results3(1,3)).*V_P_ - Q_.*trapz(t,S_t_in(t)) + Q_.*trapz(t,results3(:,3));
        %Because S_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance 
end

for i = 1:n

    g_t_in = @(t) g_in_0.*18.016 + (t > 0).*(g_up_vals(i) - g_in_0).*18.016;
    [~,results4] = simulate_alphaBetaModel_perfusion_delta(params4_,g_t_in,G_t_in,I_t_in,S_t_in,t);
    
    R_I4(i) = (results4(end,1) - results4(1,1)).*V_P_ - Q_.*trapz(t,I_t_in(t)) + Q_.*trapz(t,results4(:,1));
        %Because I_t_in is a function of time, need to calculate insulin
        %secretion by integrating the mass balance
    R_G4(i) = (results4(end,2) - results4(1,2)).*V_P_ - Q_.*trapz(t,G_t_in(t)) + Q_.*trapz(t,results4(:,2));
        %Because G_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance

    R_S4(i) = (results4(end,3) - results4(1,3)).*V_P_ - Q_.*trapz(t,S_t_in(t)) + Q_.*trapz(t,results4(:,3));
        %Because S_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance 
end
fprintf("Done")

% Insulin, Glucagon, and Somatostatin Secreted
figure(1)
subplot(3,1,1)
plot(g_up_vals,R_I, 'LineWidth',2)
title('Normal Secretion Rate of Insulin')
xlabel('Glucose (mM)')
ylabel('Insulin Secreted (mg)')
subplot(3,1,2)
plot(g_up_vals,R_G,'LineWidth',2)
title('Normal Secretion Rate of Glucagon')
xlabel('Glucose (mM)')
ylabel('Glucagon Secreted (mg)')
subplot(3,1,3)
plot(g_up_vals,R_S,'LineWidth',2)
title('Normal Secretion Rate of Somatostatin')
xlabel('Glucose (mM)')
ylabel('Somatostatin Secreted (mg)')

% glucagon secretion w multiples of somatostatin levels
% multiples of basal levels - glucagon
figure(2)
plot(g_up_vals, R_G, 'LineWidth',2)
hold on
plot(g_up_vals, R_G2, 'LineWidth',2)
plot(g_up_vals, R_G3, 'LineWidth',2)
plot(g_up_vals, R_G4, 'LineWidth',2)
hold off
% ylim([0 1.35*10^-7])
lgd = legend('Normal','2x', '3x', '0.5x');
title(lgd,'Multiples of Basal Somatostatin Levels')
title('Secretion Rate of Glucagon Varying Basal Somatostatin Levels')
xlabel('Glucose (mM)')
ylabel('Glucagon Secreted (mg)')


