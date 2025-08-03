function [t,y] = simulate_alphaBetaModel_perfusion_delta(params,g_t_in,G_t_in,I_t_in,S_t_in,t)

%simulate_alphaBetaModel_perfusion_delta runs a simulation of the complete
%alpha-cell, beta-cell, and delta-cell model in a perfusion setting. It takes in the
%parameters for the model as a vector, the glucose in flow rate trajectory
%as a function, the glucagon in flow rate trajecotry as a function, the
%insulin in flow rate trajecotry as a function, somatostatin inflow rate as a function,
% and the time values that results are desired for as a vector.
%It returns the time values that have results as a vector and the results
%of the simulation, including net signals and secretion rates that were
%calculated after the simulation, as an array.

    %Unpack parameters
    params = num2cell(params);

    [gba_, Gba_, Iba_, Sba_,...
     ~, ~, ~, ~, ~, ...
     m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
     h_IA_, n_IA_, X_A0_, m_g_, ...
     m_I_, h_I_, n_I_, ...
     m_G_, h_G_, n_G_, ...
     m_S_, h_S_,n_S_, ...
     m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
     m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
     m_S1_, h_S1_, n_S1_, m_S2_, h_S2_, n_S2_, ...
     n_SB_,h_SB_,n_SA_,h_SA_, ...
     Q_,~] = params{:};

    params = cell2mat(params);

    %Initial glucose concentrations at steady-state from mass balance
    g_0 = g_t_in(0);

    %Initial glucose signal intensities
    X_gB_0 = g_0./gba_;
    X_gA_0 = g_0./gba_;

    %Need to solve for insulin and glucose concentration so that they are
    %at steady state
    
    opts = optimoptions('fsolve','Display','none','FunctionTolerance',1e-20,'OptimalityTolerance',1e-20);
    C_0 = fsolve(@(x) SS_equations(x,params,g_t_in,G_t_in,I_t_in,S_t_in),[Gba_,Iba_,Sba_],opts);

    G_0 = C_0(1); I_0 = C_0(2); S_0 = C_0(3);

    %Initial insulin, glucagon, somatostatin signal intensities
    X_G_0 = G_0./Gba_;
    X_I_0 = I_0./Iba_;
    X_S_0 = S_0./Sba_;

    %Initial net signals
    X_B_0 = Y_B(X_gB_0,X_G_0,m_GB_,h_GB_,n_GB_,h_gB_,n_gB_,X_B0_,X_S_0,h_SA_,n_SA_);
    X_A_0 = Y_A(X_gA_0,X_I_0,h_IA_,n_IA_,X_A0_,m_g_,X_S_0,h_SA_,n_SA_);
    X_D_0 = X_I_0 + X_G_0;

    %Initial steady-state secretion rates
    R_Iss_0 = R_Iss(X_B_0,m_I_,h_I_,n_I_);
    R_Gss_0 = R_Gss(X_A_0,m_G_,h_G_,n_G_);
    R_Sss_0 = R_Sss(X_D_0,m_S_,h_S_,n_S_);

    %Initial pool masses 
    I_1_0 = R_Iss_0./hill(X_B_0,m_I1_,h_I1_,n_I1_);
    I_2_0 = hill(X_B_0,m_I1_,h_I1_,n_I1_).*I_1_0./hill(X_B_0,m_I2_,h_I2_,n_I2_);

    G_1_0 = R_Gss_0./hill(X_A_0,m_G1_,h_G1_,n_G1_);
    G_2_0 = hill(X_A_0,m_G1_,h_G1_,n_G1_).*G_1_0./hill(X_A_0,m_G2_,h_G2_,n_G2_);

    S_1_0 = R_Sss_0./hill(X_S_0,m_S1_, h_S1_, n_S1_);
    S_2_0 = hill(X_S_0,m_S1_,h_S1_,n_S1_).*S_1_0./hill(X_S_0,m_S2_,h_S2_,n_S2_);

    %Store initial condititions and time values for integration
    U0_ = [I_0;G_0;S_0; ...
           X_gB_0;X_G_0;I_1_0;I_2_0; ...
           X_gA_0;X_I_0;G_1_0;G_2_0;g_0; ...
           X_S_0; S_1_0; S_2_0;];
    t0 = min(t);
    tMax = max(t);

    %Solve system of odes and evaluate at t
    sol = ode23s(@(t,y) alphaBetaModel_perfusion_delta(t,y,params,g_t_in,G_t_in,I_t_in,S_t_in) ,[t0 tMax], U0_);
    y = deval(sol,t)';

    %Store signals and concentrations to calculate net signals and secretion
    %rates
    %Concentrations
    I = y(:,1);
    G = y(:,2);
    S = y(:,3);

    %Beta cell signals
    X_gB = y(:,4);
    X_G = y(:,5);

    %Alpha cell signals
    X_gA = y(:,8);
    X_I = y(:,9);

    % Delta cell signals
    X_S = y(:,13);

    %Calculate X_B, X_A, X_D, R_I, R_G, R_S

    %Make placeholder arrays
    X_B = zeros(length(t),1);
    X_A = zeros(length(t),1);
    X_D = zeros(length(t),1);

    for i = 1:length(t) %Cycle over time values to calculate signals and 
                        %secretion at each time
        X_B(i) = Y_B(X_gB(i),X_G(i),m_GB_,h_GB_,n_GB_,h_gB_,n_gB_,X_B0_,X_S(i),h_SB_,n_SB_);
        X_A(i) = Y_A(X_gA(i),X_I(i),h_IA_,n_IA_,X_A0_,m_g_,X_S(i),h_SA_,n_SA_);
        X_D(i) = X_G(i) + X_I(i) ; 
    end

    %This is the measured insulin/glucagon out flow rate
        %Will need to subtract out I_t_in and G_t_in to isolate the insulin
        %or glucagon secretion rate from only beta or alpha cells.
    R_I = Q_.*I;
    R_G = Q_.*G;
    R_S = Q_.*S;

    %Store the calculated net signals and secretion with the results of the
    %system of odes to return from this function
    y = [y X_B R_I X_A R_G X_D R_S];


end


% Additional functions
function s = R_Gss(X_a,m_a,h_a,n_a)
    %R_Gss represents the steady-state glucagon secretion function
    s = hill(X_a,m_a,h_a,n_a); %mg/min/islet
end 

function s = R_Iss(X_b,m_b,h_b,n_b)
    %R_Iss represents the steady-state insulin secretion function
    s = hill(X_b,m_b,h_b,n_b); %mg/min/islet
end

function s = R_Sss(X_d,m_d,h_d,n_d)
    %R_Sss represents the steady-state somatostatin secretion function
    s = hill(X_d,m_d,h_d,n_d); %mg/min/islet
end

function s = Y_A(X_gA,X_I,h_IA,n_IA,X_A0,m_g,X_S,h_SA,n_SA)
    %Y_A represents the net alpha cell signal function
    s = X_gA - hill(X_I,m_g*X_gA+X_A0,h_IA,n_IA) + X_A0 - hill(X_S,1,h_SA,n_SA);
end
    
function s = Y_B(X_gB,X_G,m_GB,h_GB,n_GB,h_gB,n_gB,X_B0,X_S,h_SB,n_SB)
    %Y_B represents the net beta cell signal function
    s = X_gB + hill(X_G,m_GB,h_GB,n_GB)*hill(X_gB,1,h_gB,n_gB) + X_B0 - hill(X_S,1,h_SB,n_SB);
end


function hi = hill(x,m,h,n)

    %Hill Function

    hi = (x >= 0) .* m./((h./x).^n + 1) + (x < 0) .* 0;
    %If the x value is less than 0, the Hill function should still be 0

end

function r = SS_equations(x,params,g_t_in,G_t_in,I_t_in,S_t_in)

    params = num2cell(params);

    [gba_, Gba_, Iba_, Sba_,...
     ~, ~, ~, ~, ...
     m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
     h_IA_, n_IA_, X_A0_, m_g_, ...
     m_I_, h_I_, n_I_, ...
     m_G_, h_G_, n_G_, ...
     m_S_, h_S_, n_S_, ...
     ~, ~, ~, ~, ~, ~, ...
     ~, ~, ~, ~, ~, ~, ...
     n_SB_, h_SB_, n_SA_, h_SA_, ...
     Q_,~] = params{:};

    %Solving for insulin and glucagon and somatostatin concentrations
    G = x(1);
    I = x(2);
    S = x(3);
    
    %Initial glucose concentrations at steady-state from mass balance
    g_0 = g_t_in(0);

    %Initial glucose signal intensities
    X_gB_0 = g_0./gba_;
    X_gA_0 = g_0./gba_;

    %Signal transduction
    X_G = G./Gba_;
    X_I = I./Iba_;
    X_S = S./Sba_;

    %Net signals
    X_B = Y_B(X_gB_0,X_G,m_GB_,h_GB_,n_GB_,h_gB_,n_gB_,X_B0_,X_S,h_SB_,n_SB_);
    X_A = Y_A(X_gA_0,X_I,h_IA_,n_IA_,X_A0_,m_g_,X_S,h_SA_,n_SA_);
    X_D = X_I + X_G;

    %Steady-state secretion - at SS, only need this, not transient portion
    R_Iss_ = R_Iss(X_B,m_I_,h_I_,n_I_);
    R_Gss_ = R_Gss(X_A,m_G_,h_G_,n_G_);
    R_Sss_ = R_Sss(X_D,m_S_,h_S_,n_S_);
    
    %Mass balances - perfusion
    r(1) = Q_.*(I_t_in(0) - I) + R_Iss_; %Perfusion 
    r(2) = Q_.*(G_t_in(0) - G) + R_Gss_; %Perfusion 
    r(3) = Q_.*(S_t_in(0) - S) + R_Sss_; %Perfusion 

end

