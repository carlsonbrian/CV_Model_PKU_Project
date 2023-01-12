% ***********************************************************************************
%         P K U / U M I C H  C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************
%
%   This function simulates the cardiovascular hemodynamics using a simple 
%   cardiovascular systems model incorporating an implementation of the Lumens
%   et al. TriSeg model.
%
%   Model originally created on     10 January 2023
%   Model last modfied on           11 January 2023
%
%   Based on the code by   Ben Randall
%                          Salla Kim
%                          Andrew Meyer
%                          Dan Beard
%
%   Heart model from       Joost Lumens et al.
%
%   Implemented by         Brian Carlson
%                          Physiological Systems Dynamics Laboratory
%                          Department of Molecular and Integrative Physiology
%                          University of Michigan
%  
% ***********************************************************************************
%  START OF  	              C V   S Y S T E M S   M O D E L   W / T R I S E G
% **********************************************************************************

    function [dxdt, outputs] = PKU_CV_dXdt(t,x,AllStruct_Struct) 


%% **********************************************************************************
%  UNPACK STRUCTS FOR         C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************

    InputData_Struct = AllStruct_Struct.InputData_Struct;
    Param_Struct = AllStruct_Struct.Param_Struct;

%% **********************************************************************************
%  UNPACK PARAMS FOR          C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************

    % Compliances (mL mmHg^(-1))
    C_SA = Param_Struct.C_SA;                       % Systemic arterial compliance
    C_SV = Param_Struct.C_SV;                       % Systemic venous compliance
    C_PA = Param_Struct.C_PA;                       % Pulmonary arterial compliance
    C_PV = Param_Struct.C_PV;                       % Pulmonary venous compliance
    % Resistances (mmHg s mL^(-1))
    R_SA  = Param_Struct.R_SA;                      % Systemic arterial resistance
    R_tSA = Param_Struct.R_tSA;                     % Trnsmrl syst art resistance
    R_PA  = Param_Struct.R_PA;                      % Pulmonary arterial resistance
    R_tPA = Param_Struct.R_tPA;                     % Trnsmrl pulm art resistance
    R_m = Param_Struct.R_m;                         % Mitral valve resistance
    R_a = Param_Struct.R_a;                         % Aortic valve resistance
    R_t = Param_Struct.R_t;                         % Tricuspid valve resistance
    R_p = Param_Struct.R_p;                         % Pulmonary valve resistance 
    % Pericardium parameters 
    Vh0 = Param_Struct.Vh0;                         % Pericrdl zero press vol (mL) 
    s   = Param_Struct.s;                           % Pericrdl stiffness coeff (uls)        
    % Force scaling factors (unitless) 
    k_pas_LV = Param_Struct.k_pas_LV;               % Scale fact passive LV force
    k_pas_RV = Param_Struct.k_pas_RV;               % Scale fact passive RV force
    k_act_LV = Param_Struct.k_act_LV;               % Scale fact active LV force
    k_act_RV = Param_Struct.k_act_RV;               % Scale fact active RV force
    % Sarcomere length parameters (um)
    Lsref   = Param_Struct.Lsref;
    Lsc0    = Param_Struct.Lsc0; 
    Lse_iso = Param_Struct.Lse_iso; 
    % Sarcomere length shortening velocity (um s^(-1))
    v_max = Param_Struct.v_max; 
    % Midwall reference surface areas (cm^2)
    Amref_LV  = Param_Struct.Amref_LV;              % Midwall reference SA for LV
    Amref_SEP = Param_Struct.Amref_SEP;             % Midwall reference SA for sept
    Amref_RV  = Param_Struct.Amref_RV;              % Midwall reference SA for RV
    % Midwall volumes (mL)
    Vw_LV  = Param_Struct.Vw_LV;                    % Midwall enclsd volume for LV
    Vw_SEP = Param_Struct.Vw_SEP;                   % Midwall enclsd volume for sept
    Vw_RV  = Param_Struct.Vw_RV;                    % Midwall enclsd volume for RV
    % Passive stress steepness parameter (dimensionless) 
    gamma = Param_Struct.gamma; 
    % Percentage of cardiac cycle (uls)
    k_TS = Param_Struct.k_TS;                       % Frac T to max systole (uls)
    k_TR = Param_Struct.k_TR;                       % Frac T max syst to basln (uls)


%% **********************************************************************************
%  UNPACK PAT DATA FOR        C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************

    HR = InputData_Struct.HR_RHCRest; 

%% **********************************************************************************
%  UNPACK VARS FOR            C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************
   
    % Axial distance of midwall junction (cm)
    xm_LV  = x(1); 
    xm_SEP = x(2); 
    xm_RV  = x(3);
    % Radial distance of midwall junction (cm)
    ym = x(4); 
    % Contractile element length (um)
    Lsc_LV  = x(5); 
    Lsc_SEP = x(6); 
    Lsc_RV  = x(7); 
    % Volumes (mL) 
    V_LV = x(8); 
    V_RV = x(9);
    V_SA = x(10); 
    V_SV = x(11);
    V_PA = x(12); 
    V_PV = x(13); 
    
%% **********************************************************************************
%  ACT FUNCTN FOR             C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************
   
    % Calculating the times to systole and relaxation for current heart rate
    T = 60/HR;                                      % Heart period (s) 
    TS_v = k_TS * T;                                % Time to max systl 
    TR_v = k_TR * T;                                % Time frm max systl to relaxtn 
    % Now computing the ventricular activation function 
    tc_v = mod(t,T);                                % Current time in the ventricles
    if tc_v >= 0 && tc_v < TS_v                     % Vntrcl contracting in systole
        Y = 0.5*(1 - cos(pi*tc_v/TS_v)); 
    elseif tc_v >= TS_v && tc_v < TR_v + TS_v       % Vntrcl relaxing in diastole
        Y = 0.5*(1 + cos(pi*(tc_v - TS_v)/TR_v)); 
    else
        Y = 0;                                      % Vntrcl relaxed and filling
    end

    
%% **********************************************************************************
%  HEART MODEL EQNS FOR       C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************
   
    % Volume of spherical cap formed by midwall surface (mL)
    Vm_LV  = -(pi / 6) * xm_LV  * (xm_LV^2  + 3 * ym^2); 
    Vm_SEP =  (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
    Vm_RV  =  (pi / 6) * xm_RV  * (xm_RV^2  + 3 * ym^2); 
    % Surface area of midwall surface (cm^2) 
    Am_LV  = pi * (xm_LV^2  + ym^2);
    Am_SEP = pi * (xm_SEP^2 + ym^2); 
    Am_RV  = pi * (xm_RV^2  + ym^2); 
    % Curvature of midwall surface (cm^(-1))
    Cm_LV  = -2 * xm_LV  / (xm_LV^2  + ym^2);
    Cm_SEP =  2 * xm_SEP / (xm_SEP^2 + ym^2);
    Cm_RV  =  2 * xm_RV  / (xm_RV^2  + ym^2);
    % Ratio of wall thickness to midwall radius of curvature (dimensionless)
    z_LV  = 3 * Cm_LV  * Vw_LV  / (2 * Am_LV); 
    z_SEP = 3 * Cm_SEP * Vw_SEP / (2 * Am_SEP); 
    z_RV  = 3 * Cm_RV  * Vw_RV  / (2 * Am_RV);
    % Myofiber strain (dimensionless)
    eps_LV  = 0.5 * log( Am_LV  / Amref_LV  ) - (1/12) * z_LV^2  - 0.019 * z_LV^4; 
    eps_SEP = 0.5 * log( Am_SEP / Amref_SEP ) - (1/12) * z_SEP^2 - 0.019 * z_SEP^4; 
    eps_RV  = 0.5 * log( Am_RV  / Amref_RV  ) - (1/12) * z_RV^2  - 0.019 * z_RV^4; 
    % Sarcomere length (um)
    Ls_LV  = Lsref * exp(eps_LV); % From geometry of whole model
    Ls_SEP = Lsref * exp(eps_SEP); 
    Ls_RV  = Lsref * exp(eps_RV); 
    % Passive stress  
    sigma_pas_LV  =  k_pas_LV * (Ls_LV - Lsc0)^gamma; 
    sigma_pas_SEP =  k_pas_LV * (Ls_SEP - Lsc0)^gamma; 
    sigma_pas_RV  =  k_pas_RV * (Ls_RV - Lsc0)^gamma; 
    % Active stress . Cell model is state variable. If geometry model
    % bigger than cell model, positive stress. 
    sigma_act_LV  = k_act_LV * Y  * (Lsc_LV  - Lsc0) * (Ls_LV  - Lsc_LV)  / Lse_iso; 
    sigma_act_SEP = k_act_LV * Y  * (Lsc_SEP - Lsc0) * (Ls_SEP - Lsc_SEP) / Lse_iso;
    sigma_act_RV  = k_act_RV * Y  * (Lsc_RV  - Lsc0) * (Ls_RV  - Lsc_RV)  / Lse_iso;
    % Total stress 
    sigma_LV  = sigma_act_LV  + sigma_pas_LV; 
    sigma_SEP = sigma_act_SEP + sigma_pas_SEP; 
    sigma_RV  = sigma_act_RV  + sigma_pas_RV; 
    % Representative midwall tension 
    Tm_LV  = (Vw_LV  * sigma_LV  / (2 * Am_LV))  * (1 + (z_LV^2)/3  + (z_LV^4)/5); 
    Tm_SEP = (Vw_SEP * sigma_SEP / (2 * Am_SEP)) * (1 + (z_SEP^2)/3 + (z_SEP^4)/5); 
    Tm_RV  = (Vw_RV  * sigma_RV  / (2 * Am_RV))  * (1 + (z_RV^2)/3  + (z_RV^4)/5);
    % Axial midwall tension component 
    Tx_LV  = -Tm_LV  * 2 * xm_LV  * ym / (xm_LV^2  + ym^2); 
    Tx_SEP =  Tm_SEP * 2 * xm_SEP * ym / (xm_SEP^2 + ym^2); 
    Tx_RV  =  Tm_RV  * 2 * xm_RV  * ym / (xm_RV^2  + ym^2); 
    % Radial midwall tension component 
    Ty_LV  = Tm_LV  * (-xm_LV^2  + ym^2) / (xm_LV^2  + ym^2); 
    Ty_SEP = Tm_SEP * (-xm_SEP^2 + ym^2) / (xm_SEP^2 + ym^2); 
    Ty_RV  = Tm_RV  * (-xm_RV^2  + ym^2) / (xm_RV^2  + ym^2);
    % Ventricular pressure 
    P_LV = -2 * Tx_LV / ym; 
    P_RV = 2 * Tx_RV / ym; 
    
%% **********************************************************************************
%  CIRC MODEL EQNS FOR        C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************
   
    % Pericardial pressure (mmHg) 
    Vh = V_LV + V_RV; 
    P_peri = exp(s*(Vh/Vh0 - 1)); 
    % Ventricular pressures (mmHg)
    P_LV = P_peri + P_LV; 
    P_RV = P_peri + P_RV; 
    % Venous pressures (mmHg)
    P_SV = V_SV / C_SV; 
    P_PV = V_PV / C_PV;
    % Valve flows (mL s^(-1))
    Q_m = max((P_PV - P_LV) / R_m, 0); 
    Q_t = max((P_SV - P_RV) / R_t, 0); 
    % When aortic valve is closed 
    Q_a = 0; 
    Q_SA = (V_SA - C_SA*P_SV)/(C_SA*(R_SA + R_tSA));
    P_SA = (R_SA*V_SA + C_SA*P_SV*R_tSA)/(C_SA*(R_SA + R_tSA));
    % When aortic valve is open 
    if (P_SA < P_LV) * (V_LV > 0)    
        Q_a = -(R_SA*V_SA - C_SA*P_LV*R_SA - C_SA*P_LV*R_tSA + ...
            C_SA*P_SV*R_tSA)/(C_SA*(R_SA*R_a + R_SA*R_tSA + R_a*R_tSA));
        Q_SA = (R_a*V_SA - C_SA*P_SV*R_a + C_SA*P_LV*R_tSA - ...
            C_SA*P_SV*R_tSA)/(C_SA*(R_SA*R_a + R_SA*R_tSA + R_a*R_tSA));
        P_SA = (R_SA*R_a*V_SA + C_SA*P_LV*R_SA*R_tSA + ...
            C_SA*P_SV*R_a*R_tSA)/(C_SA*(R_SA*R_a + R_SA*R_tSA + R_a*R_tSA));
    end
    Q_a = max(Q_a, 0); 
    % When pulmonary valve is closed 
    Q_p = 0; 
    Q_PA = (V_PA - C_PA*P_PV)/(C_PA*(R_PA + R_tPA)); 
    P_PA = (R_PA*V_PA + C_PA*P_PV*R_tPA)/(C_PA*(R_PA + R_tPA)); 
    % When pulmonary valve is open 
    if (P_PA < P_RV) * (V_RV > 0)
        Q_p = -(R_PA*V_PA - C_PA*P_RV*R_PA + C_PA*P_PV*R_tPA - ...
            C_PA*P_RV*R_tPA)/(C_PA*(R_PA*R_p + R_PA*R_tPA + R_p*R_tPA)); 
        Q_PA = (R_p*V_PA - C_PA*P_PV*R_p - C_PA*P_PV*R_tPA + ...
            C_PA*P_RV*R_tPA)/(C_PA*(R_PA*R_p + R_PA*R_tPA + R_p*R_tPA)); 
        P_PA = (R_PA*R_p*V_PA + C_PA*P_RV*R_PA*R_tPA + ...
            C_PA*P_PV*R_p*R_tPA)/(C_PA*(R_PA*R_p + R_PA*R_tPA + R_p*R_tPA));
    end 
    Q_p = max(Q_p, 0);


%% **********************************************************************************
%  ALL ODEs FOR               C V   S Y S T E M S   M O D E L   W / T R I S E G
% ***********************************************************************************
    
    % Heart model geometries which are solved implicitly and
    %  are actually not DEs. All xm and ym dimensions are found 
    %  such that all left hand side values go to zero. (1 - 4)
    dxm_LV  = -V_LV - 0.5 * Vw_LV - 0.5 * Vw_SEP + Vm_SEP - Vm_LV; 
    dxm_SEP = Tx_LV + Tx_SEP + Tx_RV;
    dxm_RV  = V_RV + 0.5 * Vw_RV + 0.5 * Vw_SEP + Vm_SEP - Vm_RV;
    dym     = Ty_LV + Ty_SEP + Ty_RV; 
    % Rate of change of the sarcomere lengths (5 - 7)
    dLsc_LV  = ((Ls_LV  - Lsc_LV)  / Lse_iso - 1) * v_max;
    dLsc_SEP = ((Ls_SEP - Lsc_SEP) / Lse_iso - 1) * v_max;
    dLsc_RV  = ((Ls_RV  - Lsc_RV)  / Lse_iso - 1) * v_max;
    % Rate of change of the circulatory compartment volumes (8 - 13)
    dV_LV = Q_m  - Q_a; 
    dV_SA = Q_a  - Q_SA; 
    dV_SV = Q_SA - Q_t; 
    dV_RV = Q_t  - Q_p; 
    dV_PA = Q_p  - Q_PA; 
    dV_PV = Q_PA - Q_m; 
    
    dxdt = [dxm_LV; dxm_SEP; dxm_RV; dym;
        dLsc_LV; dLsc_SEP; dLsc_RV; 
        dV_LV; dV_RV; dV_SA; dV_SV; dV_PA; dV_PV; 
        ]; 
    
    outputs = [P_LV; P_SA; P_SV; P_RV; P_PA; P_PV;          % 1-6
        Vm_LV; Vm_SEP; Vm_RV;                               % 7-9
        Am_LV; Am_SEP; Am_RV;                               % 10-12
        Cm_LV; Cm_SEP; Cm_RV;                               % 13-14
        eps_LV; eps_SEP; eps_RV;                            % 15-18
        sigma_pas_LV; sigma_pas_SEP; sigma_pas_RV;          % 19-21
        sigma_act_LV; sigma_act_SEP; sigma_act_RV;          % 22-24
        sigma_LV; sigma_SEP; sigma_RV;                      % 25-27
        Q_m; Q_a; Q_t; Q_p;                                 % 28-31
        Q_SA; Q_PA;                                         % 32-33
        Tm_LV; Tm_SEP; Tm_RV;                               % 34-36
        Y;                                                  % 37
        ];

end 