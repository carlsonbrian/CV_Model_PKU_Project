% ***********************************************************************************
%            P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************
%
%   This script loads data for a given patient and then simulates a cardiovacular
%   systems model which can be hand tuned to match the data
%
%   Model originally created on     10 January 2023
%   Model last modfied on           10 January 2023
%
%   Based on the code by   Ben Randall
%                          Salla Kim
%                          Andrew Meyer
%                          Dan Beard
%   Heart model from       Joost Lumens et al.
%   Implemented by         Brian Carlson
%                          Physiological Systems Dynamics Laboratory
%                          Department of Molecular and Integrative Physiology
%                          University of Michigan
%  
% ***********************************************************************************
%  START OF  	      P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************

tic

%% **********************************************************************************
%  INPUT DATA FOR  	  P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************
    
    % Patient general data
    DataSrc = "PKU";                                % Where data is from
    PID = 1;                                        % Study patient number
    Hgt = 171.8;                                    % Not given (Ave chinese male,cm)
    Wgt = 77.7;                                     % Not given (Ave chinese male,kg)
    Sex = 'M';                                      % Patient sex at birth
    % Right heart cath data
    HR_RHCRest = 54;                                % Heart rate rest RHC (beats/min)
    CO_RHCRest = 6.7;                               % Cardiac output rest RHC (L/min)
    P_SAdiast_RHCRest = 80;                         % Not given, SPbar (mmHg)
    P_SAsyst_RHCRest = 120;                         % Not given, DPbar (mmHg)
    % Transthoracic echo data
    LV_EDD_EchoRest = 48.2;                         % Use Teich Eqn, EDV_LV (mm)
    LV_ESD_EchoRest = 28.3;                         % Use Teich Eqn, ESV_LV (mm)
    V_LVdiast_EchoRest = -1;                        % Not given, MOD LV V diast (mL)
    V_LVsyst_EchoRest = -1;                         % Not given, MOD LV V syst (mL)
    % If the LV volume is not calculated with 2D Simpson's method
    %  then use Teichholz formula with the 1D top of the ventricle
    %  diameter to estimate LV volume in systole and diastole
    if (V_LVsyst_EchoRest == -1)
        % Calculate the LV systolic volume using the Teichholz expression
        V_LVsyst_EchoRest = ((LV_EDD_EchoRest/10)^3) / ...       
            ((6/pi) * ((0.075 * (LV_EDD_EchoRest/10)) + 0.18));
        V_LVdiast_EchoRest = ((LV_ESD_EchoRest/10)^3) / ...       
            ((6/pi) * ((0.075 * (LV_ESD_EchoRest/10)) + 0.18));
    end

    % Create input data structure to pass to functions
    InputData_Values = {DataSrc PID Hgt Wgt Sex HR_RHCRest ...
        CO_RHCRest P_SAdiast_RHCRest P_SAsyst_RHCRest ...
        V_LVdiast_EchoRest V_LVsyst_EchoRest};
    InputData_Fields = {'DataSrc' 'PID' 'Hgt' 'Wgt' 'Sex' ...
        'HR_RHCRest' 'CO_RHCRest' 'P_SAdiast_RHCRest' ...
        'P_SAsyst_RHCRest' 'V_LVdiast_EchoRest' ...
        'V_LVsyst_EchoRest'};
    InputData_Struct = cell2struct(InputData_Values, ...
        InputData_Fields,2);

    
%% **********************************************************************************
%  PARAMETERS FOR     P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************

    % Compliance (mL mmHg^(-1))
    C_SA = 1.875; 
    C_SV = 45.8333; 
    C_PA = 7.5; 
    C_PV = 12.5; 
    % Resistance (mmHg s mL^(-1))
    R_SA  = 1.368; 
    R_tSA = 0.08; 
    R_PA  = 0.144; 
    R_tPA = 0.02; 
    R_m = 0.02; 
    R_a = 0.0001; 
    R_t = 0.004; 
    R_p = 0.0001;  
    % Pericardium parameters 
    Vh0 = 312.5; % mL 
    s   = 10;
    % Force scaling factors (unitless) 
    k_pas_LV = 4552.2;
    k_pas_RV = 5647.3;
    k_act_LV = 3167.1; 
    k_act_RV = 7697.5;
    % Sarcomere length parameters (um)
    Lsref   = 2;
    Lsc0    = 1.51; 
    Lse_iso = 0.04; 
    % Sarcomere length shortening velocity (um s^(-1))
    v_max = 3.5; 
    % Midwall reference surface area (cm^2)
    Amref_LV  = 102.7282;                           % Midwall reference SA for LV
    Amref_SEP = 51.3641; 
    Amref_RV  = 143.2044; 
    % Midwall volume (mL)
    Vw_LV  = 82.54; 
    Vw_SEP = 41.27; 
    Vw_RV  = 25.7725; 
    % Passive stress steepness parameter (uls) 
    gamma = 7.5; 
    % Activation function parameters
    k_TS = 0.3;                                     % Frac T to max systole (uls)
    k_TR = 0.3;                                     % Frac T max syst to basln (uls)
    % Factor to adjust fraction of stressed 
    %   versus unstressed blood volume
    SVFact = 1;                                     % Factor at 1 = 30% stresse vol

    % Create parameter structure to pass to functions
    Param_Values = {C_SA C_SV C_PA C_PV R_SA  R_tSA ...
        R_PA R_tPA R_m R_a R_t R_p Vh0 s k_pas_LV ...
        k_pas_RV k_act_LV k_act_RV Lsref Lsc0 Lse_iso ...
        v_max Amref_LV  Amref_SEP Amref_RV Vw_LV ...
        Vw_SEP Vw_RV gamma k_TS k_TR SVFact};
    Param_Fields = {'C_SA' 'C_SV' 'C_PA' 'C_PV' ...
        'R_SA'  'R_tSA' 'R_PA' 'R_tPA' 'R_m' 'R_a' ...
        'R_t' 'R_p' 'Vh0' 's' 'k_pas_LV' 'k_pas_RV' ...
        'k_act_LV' 'k_act_RV' 'Lsref' 'Lsc0' 'Lse_iso' ...
        'v_max' 'Amref_LV' 'Amref_SEP' 'Amref_RV' 'Vw_LV' ...
        'Vw_SEP' 'Vw_RV' 'gamma' 'k_TS' 'k_TR' 'SVFact'};
    Param_Struct = cell2struct(Param_Values, ...
        Param_Fields,2);

    % Create one single structure of structures to pass
    AllStruct_Values = {InputData_Struct Param_Struct};
    AllStruct_Fields = {'InputData_Struct' 'Param_Struct'};
    AllStruct_Struct = cell2struct(AllStruct_Values, ...
        AllStruct_Fields,2);


%% **********************************************************************************
%  ODE SETTINGS FOR   P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************

    % ODE tolerance 
    ODE_TOL = 1e-08;


%% **********************************************************************************
%  COMPUTED VALS FOR  P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************

    T_RHCRest = 60 / HR_RHCRest;                    % T
    SV_RHC_Rest = LV_EDD_EchoRest - ...             % Do we need this?, SV
        LV_ESD_EchoRest; 
    
    % Calculate total blood volume based on height, weight and sex.
    %  This expression is from Nadler et al. Surgery 51:224,1962.
    if (Sex == 'M')
        TotBV = ((0.3669 * (Hgt/100)^3) + (0.03219 * Wgt) + 0.6041) * 1000;
    else
        TotBV = ((0.3561 * (Hgt/100)^3) + (0.03308 * Wgt) + 0.1833) * 1000;
    end
    % The original Smith model only kept track of the stressed blood volume
    %  and ignores the unstressed volume - in essence it is assuming the unstressed
    %  volume will only be recruited due to a change in compliance in the
    %  compartments. Assuming they were simulating a typical 5000 mL total blood
    %  volume they included only 1500 mL (or 30%) in the circulating volume
    %  therefore we will multiply our calculated TotBV value by 0.3 to yield 
    %  circulating blood volume. To account for extra recruited volume in heart
    %  disease the 30% circulating blood volume can be altered by changing SVFact
    StressBV = SVFact * 0.30 * TotBV;
    
    % Calculate number of beats/time span

    
    
    % Times (s) and other holdovers from data struct
    dt = 1e-03;
    tspan = [0:dt:20];  
    T     = 60 / HR_RHCRest; 
    SPbar = 120;
    DPbar = 80;
%     Vtot = 5000;
    CO_lv = (CO_RHCRest * 1000) / 60;               % mL/s
    SV    = CO_lv*T; 
    EDV_LV = 125;
    ESV_LV = EDV_LV - SV;

%     % Unstressed volumes (mL) 
%     V_SAu = 525;
%     V_SVu = 2475; 
%     V_PAu = 100; 
%     V_PVu = 900;


%% **********************************************************************************
%  INIT CONDS FOR     P K U / U M I C H   H F p E F   P R O J E C T   S C R I P T
% ***********************************************************************************

    % Setting the initial conditions for the cardiac geometry
    xm_LV_0 = 4.6147;
    xm_SEP_0 = 0.6066;
    xm_RV_0 = 4.6698;
    ym_0 = 4.0621;
    % Setting the initial conditions for the sarcomere lengths
    Lsc_LV_0 = 2.0985;
    Lsc_SEP_0 = 1.9903;
    Lsc_RV_0 = 1.7925;
    % Setting state variable initial conditions 
    %  Note: the initial volume division is scaled as a fraction
    %  of stressed volume according to Beneken Table 1-1 where 
    %  total blood volume is 4544 mL and stressed/total is 18.75%
    V_LV_0 = (175/852) * StressBV;                  % Dan/Andrew value 125 mL
    V_RV_0 = (175/852) * StressBV;                  % Dan/Andrew value 125 mL
    V_SA_0 = (160/852) * StressBV;                  % Dan/Andrew value 150 mL
    V_SV_0 = (219/852) * StressBV;                  % Dan/Andrew value 45.83 mL
    V_PA_0 = (69/852)  * StressBV;                  % Dan/Andrew value 75 mL
    V_PV_0 = (54/852)  * StressBV;                  % Dan/Andrew value 25 mL

    X0 = [xm_LV_0; xm_SEP_0; xm_RV_0; ym_0;                 % 1-4
        Lsc_LV_0; Lsc_SEP_0; Lsc_RV_0;                      % 5-7
        V_LV_0; V_RV_0; V_SA_0; V_SV_0; V_PA_0; V_PV_0;     % 8-13
        ]; 


    %% Set mass matrix M for DAE 

    M = speye(length(init));
    M(1,1) = 0;
    M(2,2) = 0;
    M(3,3) = 0;
    M(4,4) = 0; 

    %% Solve model 
    % Use a while loop to allow model to converge to steady-state 
    
    ndone = 0; 
    while ndone == 0 
        
        % Set options and compute model for time span tspan 
        opts = odeset('Mass',M,'RelTol',ODE_TOL,'AbsTol',ODE_TOL);
        sol  = ode15s(@PKU_CV_dXdt,tspan,X0,opts,AllStruct_Struct);
       
    
        if sol.x(end) ~= tspan(end) 
            % Check to see if the model solved to the end of tspan. If not,
            % set the initial conditions for the next loop at the start of 
            % the previous beat and solve for 10 more beats 

            t = sol.x(1):dt:sol.x(end); 
            beats = mod(t,T); 
            x = find(round(beats,3) == 0);
            y = find(t(x) <= t(end)); 
            tspan = tspan(1):dt:tspan(x(y(end)));
        
            sols = deval(sol,tspan);
           
            if round(x(y(end-1))) == 1 
                % If model stopped after one cardiac cycle, reset initial
                % conditions to the beginning of the second cardiac cycle.
                init  = sols(:,x(y(end))); 
                tspan = tspan(x(y(end))):dt:tspan(x(y(end))) + 10*T;
            else 
                % If the model stopped after more than one cardiac cycle, 
                % reset the initial conditions to the previous fully 
                % computed cardiac cycle. 
                init  = sols(:,x(y(end-1))); 
                tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
            end 
             
        else 
            % If the model has successfully solved at least 10 beats, then we can
            % assess whether the model has reached steady state 
            
            % Extract sarcomere lengths and systemic arterial pressure (P_SA)     
            sols = deval(sol,tspan);
                
            Lsc_LV  = sols(5,:);
            Lsc_SEP = sols(6,:); 
            Lsc_RV  = sols(7,:); 
            P_SA    = sols(10,:) / C_SA; 
            

            
            % Find the last 5 beats of the simulation 
            xx = find(tspan >= tspan(end) - 5*T); 
            
            % Set a peak threshold as half of the amplitude of the last 5 beats 
            threshold_LV  = (max(Lsc_LV(xx))  - min(Lsc_LV(xx)))/2;
            threshold_SEP = (max(Lsc_SEP(xx)) - min(Lsc_SEP(xx)))/2;
            threshold_RV  = (max(Lsc_RV(xx))  - min(Lsc_RV(xx)))/2;
            
            % Determine the length of half of the cardiac cycle 
            half_per = round((T/2)/dt); 
            
            % Find peaks for the last 5 beats with conditions that the
            % amplitude must be at least half of the amplitude of the last
            % 5 peats and the minimum distance between peaks must be at
            % least half of the cardiac cycle apart 
            [pks_Lsc_LV, loc_pks_Lsc_LV]  = findpeaks(...
                Lsc_LV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_LV); 
            [pks_Lsc_SEP,loc_pks_Lsc_SEP] = findpeaks(...
                Lsc_SEP,'MinPeakDistance',half_per,'MinPeakProminence',threshold_SEP); 
            [pks_Lsc_RV, loc_pks_Lsc_RV]  = findpeaks(...
                Lsc_RV,'MinPeakDistance',half_per,'MinPeakProminence',threshold_RV); 
            [pks_P_SA,   loc_pks_P_SA]    = findpeaks(...
                P_SA,'MinPeakDistance',half_per); 
            
            % Exclude the last peak, so there are 4 peaks (This is
            % important in case something unexpected happens at the end of
            % the signal) 
            pks_Lsc_LV  = pks_Lsc_LV(end-5:end-1); 
            pks_Lsc_SEP = pks_Lsc_SEP(end-5:end-1); 
            pks_Lsc_RV  = pks_Lsc_RV(end-5:end-1);
            pks_P_SA    = pks_P_SA(end-5:end-1); 
            
            % Find the locations of the peaks 
            loc_pks_Lsc_LV  = loc_pks_Lsc_LV(end-5:end-1); 
            loc_pks_Lsc_SEP = loc_pks_Lsc_SEP(end-5:end-1); 
            loc_pks_Lsc_RV  = loc_pks_Lsc_RV(end-5:end-1); 
            loc_pks_P_SA    = loc_pks_P_SA(end-5:end-1); 
            
            % Find the times where the peaks occur 
            t_pks_Lsc_LV  = tspan(loc_pks_Lsc_LV);
            t_pks_Lsc_SEP = tspan(loc_pks_Lsc_SEP);
            t_pks_Lsc_RV  = tspan(loc_pks_Lsc_RV);
            t_pks_P_SA    = tspan(loc_pks_P_SA); 
            
            % Create a linear regression through the peaks 
            pf_Lsc_LV  = polyfit(t_pks_Lsc_LV,pks_Lsc_LV,1); 
            pf_Lsc_SEP = polyfit(t_pks_Lsc_SEP,pks_Lsc_SEP,1); 
            pf_Lsc_RV  = polyfit(t_pks_Lsc_RV,pks_Lsc_RV,1); 
            pf_P_SA    = polyfit(t_pks_P_SA,pks_P_SA,1); 
            
            % Extract the slope of the regression line 
            slope_Lsc_LV  = pf_Lsc_LV(1);
            slope_Lsc_SEP = pf_Lsc_SEP(1);
            slope_Lsc_RV  = pf_Lsc_RV(1);
            slope_P_SA    = pf_P_SA(1);
            
            
            % If the slope is sufficiently small (i.e. the line is 
            % practically horizontal), we have reached steady state 
            slope_lim = 1e-3;
                % Stopping criteria 
                if abs(slope_P_SA) < slope_lim && abs(slope_Lsc_LV) < slope_lim && ...
                        abs(slope_Lsc_SEP) < slope_lim && abs(slope_Lsc_RV) < slope_lim
                    ndone = 1; 
                end 
                
            % If we have not reached steady-state, solve the model for 10 more
            % beats and reassess convergence 
            beats = mod(tspan,T); 
            x = find(round(beats,3) == 0);
            y = find(tspan(x) <= tspan(end));
            tspan = tspan(x(y(end-1))):dt:tspan(x(y(end-1))) + 10*T;
            init  = sols(:,x(y(end-1))); 
        end 
    end
    
    % After determining that the model is in steady-state, solve 2 more beats 
    time = 0:dt:2*T; 
    sol  = ode15s(@PKU_CV_dXdt,[time(1) time(end)],init,opts,AllStruct_Struct);
    sols = deval(sol,time);
    sols = sols'; 
    
    %% Calculate other time-varying model quantities using auxiliary equations
    
    o = zeros(37,length(time));  
    for i = 1:length(time) 
        [~,o(:,i)] = PKU_CV_dXdt(time(i),sols(i,:),AllStruct_Struct);
    end 
    
    %% Outputs 
    
    %outputs.time = time; 
    
    xm_LV  = sols(:,1); 
    xm_SEP = sols(:,2); 
    xm_RV  = sols(:,3); 
    ym     = sols(:,4); 
    
    Lsc_LV  = sols(:,5);
    Lsc_SEP = sols(:,6); 
    Lsc_RV  = sols(:,7); 
    
    V_LV = sols(:,8); 
    V_RV = sols(:,9); 
    V_SA = sols(:,10) + V_SAu; 
    V_SV = sols(:,11) + V_SVu; 
    V_PA = sols(:,12) + V_PAu; 
    V_PV = sols(:,13) + V_PVu; 
    Vtot = sum(sols(end,8:13)) + V_SAu + V_SVu + V_PAu + V_PVu; 
    
    P_LV = o(1,:)'; 
    P_SA = o(2,:)'; 
    P_SV = o(3,:)'; 
    P_RV = o(4,:)'; 
    P_PA = o(5,:)'; 
    P_PV = o(6,:)'; 
    
    Vm_LV  = o(7,:)'; 
    Vm_SEP = o(8,:)'; 
    Vm_RV  = o(9,:)'; 
    
    Am_LV  = o(19,:)'; 
    Am_SEP = o(11,:)'; 
    Am_RV  = o(12,:)'; 
    
    Cm_LV  = o(13,:)';
    Cm_SEP = o(14,:)';
    Cm_RV  = o(15,:)'; 
    
    eps_LV  = o(16,:)'; 
    eps_SEP = o(17,:)'; 
    eps_RV  = o(18,:)'; 
    
    sigma_pas_LV  = o(19,:)';
    sigma_pas_SEP = o(20,:)';
    stresses.passive.sigma_pas_RV  = o(21,:)';
    
    sigma_act_LV  = o(22,:)';
    sigma_act_SEP = o(23,:)';
    sigma_act_RV  = o(24,:)';
    % Totals
    sigma_LV  = o(25,:)';
    sigma_SEP = o(26,:)';
    sigma_RV  = o(27,:)';
    
    % Convert m^3 s^(-1) to L min^(-1)
    Q_m = o(28,:)' * 1e-3 * 60; 
    Q_a = o(29,:)' * 1e-3 * 60; 
    Q_t = o(30,:)' * 1e-3 * 60; 
    Q_p = o(31,:)' * 1e-3 * 60; 
    
    Q_SA = o(32,:)' * 1e-3 * 60; 
    Q_PA = o(33,:)' * 1e-3 * 60; 
   
    Tm_LV  = o(34,:)';
    Tm_SEP = o(35,:)'; 
    Tm_RV  = o(36,:)'; 

    Y = o(37,:)'; 
    
    %% Find times for end-systole and end-diastole 

    i_ES = find(diff(Q_m) > 0,1,'first'); 
    i_ED = find(diff(Q_a) > 0,1,'first'); 
    
    ES = time(i_ES); 
    ED = time(i_ED); 

    %% Figures 

    purple = [148, 0, 211]/255; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PV loop 
    hfig1 = figure(1); 
    clf
    hold on 
    h1 = plot(V_LV, P_LV, 'r','linewidth',2);
    h2 = plot(V_RV, P_RV, 'b','linewidth',2);
    xlabel('Volume (mL)')
    ylabel('Pressure (mmHg)')
    legend([h1 h2],'LV','RV')

    x_max = max(vertcat(V_LV,V_RV));
    x_max = max(x_max,150); 
    y_max = max(vertcat(P_LV,P_RV)); 
    xlim([0,round(x_max)+10])
    ylim([0,round(y_max)+10])
    set(gca,'FontSize',20)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Psa vs time 
    hfig2 = figure(2);  
    clf
    hold on 
    plot([time(1) time(end)],SPbar * ones(2,1),'k:','linewidth',0.5)
    plot([time(1) time(end)],DPbar * ones(2,1),'k:','linewidth',0.5)
    h1 = plot(time,P_SA, 'color', purple, 'linewidth',2); 
    xlabel('Time (s)')
    ylabel('Pressure (mmHg)','linewidth',2)
    legend([h1],'P_{SA}')
    ylim([round(DPbar)-20 round(SPbar)+20])
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % High pressures vs time 
    hfig3 = figure(3);
    clf 
    hold on
    yline(SPbar,':')
    yline(DPbar,':')
    h1 = plot(time,P_LV,'r','linewidth',2);
    h2 = plot(time,P_SA,'color',purple,'linewidth',2);
    legend([h1 h2],'P_{LV}','P_{SA}','orientation','horizontal')
    xlabel('Time (s)')
    ylabel('Pressure (mmHg)')
    
    y_max = max(vertcat(P_LV,P_SA)); 
    xlim([time(1) time(end)])
    ylim([0 round(y_max)+10])
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Low pressures vs time 
    hfig4 = figure(4);
    clf 
    hold on
    plot(time,P_SV,'k','linewidth',2)
    plot(time,P_RV,'b','linewidth',2)
    plot(time,P_PA,'c','linewidth',2)
    plot(time,P_PV,'m','linewidth',2)
    legend('P_{SV}','P_{RV}','P_{PA}','P_{PV}','orientation','horizontal')
    xlabel('Time (s)')
    ylabel('Pressure (mmHg)')
    set(gca,'FontSize',20)
    xlim([time(1) time(end)])
    ylim([0 max(P_PA + 10)])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ventricular volumes vs time 
    hfig5 = figure(5); 
    clf
    hold on 
    yline(EDV_LV,':')
    yline(ESV_LV,':')
    h1 = plot(time,V_LV,'r','linewidth',2);
    h2 = plot(time,V_RV,'b','linewidth',2);
    legend([h1, h2],'LV','RV','orientation','horizontal')
    xlabel('Time (s)')
    ylabel('Volume (mmHg)')
    set(gca,'FontSize',20)
    xlim([time(1) time(end)])
    if max(V_LV) >= max(V_RV)
        maxV = max(V_LV);
    else 
        maxV = max(V_RV);
    end
    if min(V_LV) <= min(V_RV)
        minV = min(V_LV);
    else 
        minV = min(V_RV);
    end
    ylim([(minV-50) (maxV+50)])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Valve flows vs time 
    hfig6 = figure(6); 
    clf
    hold on 
    plot(time,Q_t,'b','linewidth',2)
    plot(time,Q_p,'b--','linewidth',2)
    plot(time,Q_a,'r--','linewidth',2)
    plot(time,Q_m,'r','linewidth',2)
    legend('Q_{t}','Q_{p}','Q_{a}','Q_{m}')
    xlabel('Time (s)')
    ylabel('Flow (L min^{-1})')
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sarcomere lenghts vs time 
    hfig7 = figure(7); 
    clf
    hold on 
    plot(time,Lsc_LV,'r')
    plot(time,Lsc_SEP,'g')
    plot(time,Lsc_RV,'b')
    legend('LV','SEP','RV')
    xlabel('Time (s)')
    ylabel('Sarcomere length (\mum)')
    set(gca,'FontSize',20)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activation function 
    hfig8 = figure(8); 
    clf
    hold on 
    xline(ES,':')
    xline(ES+T,':')
    xline(ED,'--')
    xline(ED+T,'--')
    plot(time,Y,'k','linewidth',3)
    set(gca,'FontSize',20)
    xlabel('Time (s)')
    ylabel('Y(t)')
    title('Cardiac Activation')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Septal curvature 
    hfig9 = figure(9); 
    clf
    hold on
    xline(ES,':')
    xline(ES+T,':')
    xline(ED,'--')
    xline(ED+T,'--')
    plot(time,Cm_SEP,'k','linewidth',2) 
    ylabel('Septal Curvature (cm^{-1})')
    xlabel('Time (s)')
    set(gca,'FontSize',20)
    ylim([0 .5])

    toc



