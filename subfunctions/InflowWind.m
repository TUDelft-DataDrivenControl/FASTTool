function InflowWind(Wind,U,H,R,windSeed)

% Wind types:
% 1 Steady wind
% 2 Stepped wind
% 3 Normal wind profile (NWP)
% 4 Normal turbulence model (NTM)
% 5 Annual extreme wind speed (EWM)
% 6 50-year extreme wind speed (EWM)
% 7 Extreme wind shear (EWS)
% 8 Extreme turbulence model (ETM)
% 9 Extreme operating gust (EOG)
% 10 Extreme direction change (EDC)
% 11 Extreme coherent gust (ECG)

%% InflowWind input file including steady wind (1) and NWP (3)
if Wind.Type == 1
    WindType = 1;
    PLexp = 0;
elseif Wind.Type == 3
    WindType = 1;
    PLexp = 0.2;
else
    WindType = 4;
    PLexp = 0.2;
end

fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'InflowWind.dat'], 'wt');
fprintf(fid, '------- InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------\n');
fprintf(fid, '12 m/s turbulent winds on 31x31 FF grid and tower for FAST CertTests #18, #19, #21, #22, #23, and #24\n');
fprintf(fid, '---------------------------------------------------------------------------------------------------------------\n');
fprintf(fid, 'False         Echo           - Echo input data to <RootName>.ech (flag)\n');
fprintf(fid, '          %i  WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined)\n', WindType);
fprintf(fid, '          0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)\n');
fprintf(fid, '          1   NWindVel       - Number of points to output the wind velocity    (0 to 9)\n');
fprintf(fid, '          0   WindVxiList    - List of coordinates in the inertial X direction (m)\n');
fprintf(fid, '          0   WindVyiList    - List of coordinates in the inertial Y direction (m)\n');
fprintf(fid, '      %5.4f   WindVziList    - List of coordinates in the inertial Z direction (m)\n', H);
fprintf(fid, '================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================\n');
fprintf(fid, '      %5.4f   HWindSpeed     - Horizontal windspeed                            (m/s)\n', U);
fprintf(fid, '      %5.4f   RefHt          - Reference height for horizontal wind speed      (m)\n', H);
fprintf(fid, '      %5.4f   PLexp          - Power law exponent                              (-)\n', PLexp);
fprintf(fid, '================== Parameters for Uniform wind file   [used only for WindType = 2] ============================\n');
fprintf(fid, '"unused"    Filename       - Filename of time series data for uniform wind field.      (-)\n');
fprintf(fid, '         90   RefHt          - Reference height for horizontal wind speed                (m)\n');
fprintf(fid, '     125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)\n');
fprintf(fid, '================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============\n');
fprintf(fid, '"unused"    Filename       - Name of the Full field wind file to use (.bts)\n');
fprintf(fid, '================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========\n');
fprintf(fid, '"wind"    FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)\n');
fprintf(fid, 'False         TowerFile      - Have tower file (.twr) (flag)\n');
fprintf(fid, '================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================\n');
fprintf(fid, '"unused"    FileName_u     - name of the file containing the u-component fluctuating wind (.bin)\n');
fprintf(fid, '"unused"    FileName_v     - name of the file containing the v-component fluctuating wind (.bin)\n');
fprintf(fid, '"unused"    FileName_w     - name of the file containing the w-component fluctuating wind (.bin)\n');
fprintf(fid, '         64   nx             - number of grids in the x direction (in the 3 files above) (-)\n');
fprintf(fid, '         32   ny             - number of grids in the y direction (in the 3 files above) (-)\n');
fprintf(fid, '         32   nz             - number of grids in the z direction (in the 3 files above) (-)\n');
fprintf(fid, '         16   dx             - distance (in meters) between points in the x direction    (m)\n');
fprintf(fid, '          3   dy             - distance (in meters) between points in the y direction    (m)\n');
fprintf(fid, '          3   dz             - distance (in meters) between points in the z direction    (m)\n');
fprintf(fid, '         90   RefHt          - reference height; the height (in meters) of the vertical center of the grid (m)\n');
fprintf(fid, '  -------------   Scaling parameters for turbulence   ---------------------------------------------------------\n');
fprintf(fid, '          1   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]\n');
fprintf(fid, '          1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]\n');
fprintf(fid, '          1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]\n');
fprintf(fid, '          1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]\n');
fprintf(fid, '         12   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]\n');
fprintf(fid, '          8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]\n');
fprintf(fid, '          2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]\n');
fprintf(fid, '  -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------\n');
fprintf(fid, '          5   URef           - Mean u-component wind speed at the reference height (m/s)\n');
fprintf(fid, '          2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)\n');
fprintf(fid, '        0.2   PLExp          - Power law exponent (-) (used for PL wind profile type only)\n');
fprintf(fid, '       0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)\n');
fprintf(fid, '====================== OUTPUT ==================================================\n');
fprintf(fid, 'False         SumPrint     - Print summary data to <RootName>.IfW.sum (flag)\n');
fprintf(fid, '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n');
fprintf(fid, '"Wind1VelX"               X-direction wind velocity at point WindList(1)\n');
fprintf(fid, '"Wind1VelY"               Y-direction wind velocity at point WindList(1)\n');
fprintf(fid, '"Wind1VelZ"               Z-direction wind velocity at point WindList(1)\n');
fprintf(fid, 'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n');
fprintf(fid, '---------------------------------------------------------------------------------------\n');
fclose(fid);

%% Stepped wind profile (2)
if Wind.Type == 2
    
    % Wind profile
    t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
    if rem(length(t),2) ~= 0
        t = [t, max(t)+Wind.dt];
    end
    x = U*t;
    y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
    z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
    u = Wind.Step-1+ceil((U-Wind.Step+1)*t/Wind.T);
    u(u<0) = 0;
    u(u>U) = U;
    u = repmat(u(:),[1,Wind.Ny,Wind.Nz]);
    v = zeros(size(u));
    w = zeros(size(u));
    
    % Save .wnd file
    writebladed([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'wind'],(u-U)/U,v/U,w/U,x,y,z,U);
    
    % Save .sum file
    fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'wind.sum'], 'wt');
    fprintf(fid, 'T\tCLOCKWISE\n');
    fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', H);
    fprintf(fid, '%0.3f\tUBAR\n', U);
    fprintf(fid, '%0.3f\tTI(u)\n', 100);
    fprintf(fid, '%0.3f\tTI(v)\n', 100);
    fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
    fprintf(fid, '0\tHEIGHT OFFSET');
    fclose(fid);
    
end

%% TurbSim for NTM (4), EWM1 (5), EWM50 (6), ETM (8)
if Wind.Type == 4 || Wind.Type == 5 || Wind.Type == 6 || Wind.Type == 8
    
    % IEC class
    ABC = 'ABC';
    IECturbc = ABC(Wind.Class(2));
    if Wind.Type == 4
        IEC_WindType = 'NTM';
    elseif Wind.Type == 5
        IEC_WindType = [int2str(Wind.Class(1)), 'EWM1'];
    elseif Wind.Type == 6
        IEC_WindType = [int2str(Wind.Class(1)), 'EWM50'];
    elseif Wind.Type == 8
        IEC_WindType = [int2str(Wind.Class(1)), 'ETM'];
    end
    
    % Turbulence seed
    if windSeed < 0
        rng('shuffle')
        seed = randi([-2147483648, 2147483647]);
    else
        seed = windSeed;
    end
    
    % Write TurbSim input file
    fid = fopen('subfunctions\inputfiles\wind.inp', 'wt');
    fprintf(fid, 'TurbSim Input File. Valid for TurbSim v1.06.00, 21-Sep-2012\n');
    fprintf(fid, '\n');
    fprintf(fid, '---------Runtime Options-----------------------------------\n');
    fprintf(fid, '%i                  RandSeed1       - First random seed  (-2147483648 to 2147483647) \n', seed);
    fprintf(fid, 'RANLUX              RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n');
    fprintf(fid, 'False               WrBHHTP         - Output hub-height turbulence parameters in binary form?  (Generates RootName.bin)\n');
    fprintf(fid, 'False               WrFHHTP         - Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat)\n');
    fprintf(fid, 'False               WrADHH          - Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh)\n');
    fprintf(fid, 'False               WrADFF          - Output full-field time-series data in TurbSim/AeroDyn form? (Generates Rootname.bts)\n');
    fprintf(fid, 'True                WrBLFF          - Output full-field time-series data in BLADED/AeroDyn form?  (Generates RootName.wnd)\n');
    fprintf(fid, 'False               WrADTWR         - Output tower time-series data? (Generates RootName.twr)\n');
    fprintf(fid, 'False               WrFMTFF         - Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, RootName.v, RootName.w)\n');
    fprintf(fid, 'True                WrACT           - Output coherent turbulence time steps in AeroDyn form? (Generates RootName.cts)\n');
    fprintf(fid, 'True                Clockwise       - Clockwise rotation looking downwind? (used only for full-field binary files - not necessary for AeroDyn)\n');
    fprintf(fid, '0                   ScaleIEC        - Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale uniformly; 2=use individual scales]\n');
    fprintf(fid, ' \n');
    fprintf(fid, '--------Turbine/Model Specifications-----------------------\n');
    fprintf(fid, '%i                  NumGrid_Z       - Vertical grid-point matrix dimension\n', Wind.Nz);
    fprintf(fid, '%i                  NumGrid_Y       - Horizontal grid-point matrix dimension\n', Wind.Ny);
    fprintf(fid, '%2.5f               TimeStep        - Time step [seconds]\n', Wind.dt);
    fprintf(fid, '%2.3f               AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n', Wind.T);
    fprintf(fid, '%2.3f               UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds)\n', Wind.T);
    fprintf(fid, '%2.3f               HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n', H);
    fprintf(fid, '%2.3f               GridHeight      - Grid height [m] \n', Wind.Lz);
    fprintf(fid, '%2.3f               GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n', Wind.Ly);
    fprintf(fid, '0                   VFlowAng        - Vertical mean flow (uptilt) angle [degrees]\n');
    fprintf(fid, '0                   HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n');
    fprintf(fid, '  \n');
    fprintf(fid, '--------Meteorological Boundary Conditions-------------------\n');
    fprintf(fid, '"IECKAI"            TurbModel       - Turbulence model ("IECKAI"=Kaimal, "IECVKM"=von Karman, "GP_LLJ", "NWTCUP", "SMOOTH", "WF_UPW", "WF_07D", "WF_14D", "TIDAL", or "NONE")\n');
    fprintf(fid, '"1-ED3"             IECstandard     - Number of IEC 61400-x standard (x=1,2, or 3 with optional 61400-1 edition number (i.e. "1-Ed2") )\n');
    fprintf(fid, '"%s"                IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP model, not used for other models)\n', IECturbc);
    fprintf(fid, '"%s"                IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n', IEC_WindType);
    fprintf(fid, 'default             ETMc            - IEC Extreme Turbulence Model "c" parameter [m/s]\n');
    fprintf(fid, '"PL"                WindProfileType - Wind profile type ("JET";"LOG"=logarithmic;"PL"=power law;"H2L"=Log law for TIDAL spectral model;"IEC"=PL on rotor disk, LOG elsewhere; or "default")\n');
    fprintf(fid, '%2.3f               RefHt           - Height of the reference wind speed [m]\n', H);
    fprintf(fid, '%2.3f               URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET wind profile)\n', U);
    fprintf(fid, 'default             ZJetMax         - Jet height [m] (used only for JET wind profile, valid 70-490 m)\n');
    fprintf(fid, 'default             PLExp           - Power law exponent [-] (or "default")           \n');
    fprintf(fid, 'default             Z0              - Surface roughness length [m] (or "default")\n');
    fprintf(fid, '\n');
    fprintf(fid, '--------Non-IEC Meteorological Boundary Conditions------------\n');
    fprintf(fid, 'default             Latitude        - Site latitude [degrees] (or "default")\n');
    fprintf(fid, '0.05                RICH_NO         - Gradient Richardson number \n');
    fprintf(fid, 'default             UStar           - Friction or shear velocity [m/s] (or "default")\n');
    fprintf(fid, 'default             ZI              - Mixing layer depth [m] (or "default")\n');
    fprintf(fid, 'default             PC_UW           - Hub mean u''w'' Reynolds stress (or "default")\n');
    fprintf(fid, 'default             PC_UV           - Hub mean u''v'' Reynolds stress (or "default")\n');
    fprintf(fid, 'default             PC_VW           - Hub mean v''w'' Reynolds stress (or "default")\n');
    fprintf(fid, 'default             IncDec1         - u-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")\n');
    fprintf(fid, 'default             IncDec2         - v-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")\n');
    fprintf(fid, 'default             IncDec3         - w-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")\n');
    fprintf(fid, 'default             CohExp          - Coherence exponent (or "default")\n');
    fprintf(fid, '\n');
    fprintf(fid, '--------Coherent Turbulence Scaling Parameters-------------------\n');
    fprintf(fid, '"dummy"             CTEventPath     - Name of the path where event data files are located\n');
    fprintf(fid, '"Random"            CTEventFile     - Type of event files ("LES", "DNS", or "RANDOM")\n');
    fprintf(fid, 'true                Randomize       - Randomize the disturbance scale and locations? (true/false)\n');
    fprintf(fid, '1.0                 DistScl         - Disturbance scale (ratio of wave height to rotor disk). (Ignored when Randomize = true.)\n');
    fprintf(fid, '0.5                 CTLy            - Fractional location of tower centerline from right (looking downwind) to left side of the dataset. (Ignored when Randomize = true.)\n');
    fprintf(fid, '0.5                 CTLz            - Fractional location of hub height from the bottom of the dataset. (Ignored when Randomize = true.)\n');
    fprintf(fid, '30.0                CTStartTime     - Minimum start time for coherent structures in RootName.cts [seconds]\n');
    fprintf(fid, '\n');
    fprintf(fid, '==================================================\n');
    fprintf(fid, 'NOTE: Do not add or remove any lines in this file!\n');
    fprintf(fid, '==================================================\n');
    fclose(fid);
    
    % Call to TurbSim
%     system('subfunctions\TurbSim [/h] "subfunctions\inputfiles\wind.inp"'); % Console output
    [~, ~] = system('subfunctions\TurbSim [/h] "subfunctions\inputfiles\wind.inp"'); % No console output
    
end

%% EWS (7)
if Wind.Type == 7
    
    % IEC class
    if Wind.Class(1) == 1
        Vref = 50;
    elseif Wind.Class(1) == 2
        Vref = 42.5;
    elseif Wind.Class(1) == 3
        Vref = 37.5;
    end
    if Wind.Class(2) == 1
        Iref = 0.16;
    elseif Wind.Class(2) == 2
        Iref = 0.14;
    elseif Wind.Class(2) == 3
        Iref = 0.12;
    end
    if H >= 60
        Lambda1 = 0.7*H;
    else
        Lambda1 = 42;
    end
    
    % NTM values
    sigma1 = Iref*(0.75*U+5.6);
    alpha = 0.2;
    
    % EWS parameters
    beta = 6.4;
    Tg = 12;
    T0 = Wind.EWS;

    % Wind profile
    t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
    if rem(length(t),2) ~= 0
        t = [t, t+Wind.dt];
    end
    x = U*t;
    y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
    z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
    [T,~,Z] = ndgrid(t,y,z);
    u = U*(Z/H).^alpha + (Z-H)/(2*R) * (2.5 + 0.2*beta*sigma1*(2*R/Lambda1).^0.25) .* (1 -cos(2*pi*(T-T0)/Tg));
    u(T<T0) = U*(Z(T<T0)/H).^alpha;
    u(T>T0+Tg) = U*(Z(T>T0+Tg)/H).^alpha;
    v = zeros(size(u));
    w = zeros(size(u));
    
    % Save .wnd file
    writebladed('subfunctions\inputfiles\wind',(u-U)/U,v/U,w/U,x,y,z,U);
    
    % Save .sum file
    fid = fopen('subfunctions\inputfiles\wind.sum', 'wt');
    fprintf(fid, 'T\tCLOCKWISE\n');
    fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', H);
    fprintf(fid, '%0.3f\tUBAR\n', U);
    fprintf(fid, '%0.3f\tTI(u)\n', 100);
    fprintf(fid, '%0.3f\tTI(v)\n', 100);
    fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
    fprintf(fid, '0\tHEIGHT OFFSET');
    fclose(fid);
    
end

%% EOG (9)
if Wind.Type == 9
    
    % IEC class
    if Wind.Class(1) == 1
        Vref = 50;
    elseif Wind.Class(1) == 2
        Vref = 42.5;
    elseif Wind.Class(1) == 3
        Vref = 37.5;
    end
    if Wind.Class(2) == 1
        Iref = 0.16;
    elseif Wind.Class(2) == 2
        Iref = 0.14;
    elseif Wind.Class(2) == 3
        Iref = 0.12;
    end
    if H >= 60
        Lambda1 = 0.7*H;
    else
        Lambda1 = 42;
    end
    
    % EWM values
    Ve50 = 1.4*Vref;
    Ve1 = 0.8*Ve50;
    
    % NTM values
    alpha = 0.2;
    sigma1 = Iref*(0.75*U+5.6);
    
    % EOG parameters
    Vgust = min([1.35*(Ve1-U), 3.3*sigma1/(1+0.1*2*R/Lambda1)]);
    Tg = 10.5;
    T0 = Wind.EOG;

    % Wind profile
    t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
    if rem(length(t),2) ~= 0
        t = [t, t+Wind.dt];
    end
    x = U*t;
    y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
    z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
    [T,~,Z] = ndgrid(t,y,z);
    u = U*(Z/H).^alpha - 0.37*Vgust*sin(3*pi*(T-T0)/Tg) .* (1-cos(2*pi*(T-T0)/Tg));
    u(T<T0) = U*(Z(T<T0)/H).^alpha;
    u(T>T0+Tg) = U*(Z(T>T0+Tg)/H).^alpha;
    v = zeros(size(u));
    w = zeros(size(u));
    
    % Save .wnd file
    writebladed('subfunctions\inputfiles\wind',(u-U)/U,v/U,w/U,x,y,z,U);
    
    % Save .sum file
    fid = fopen('subfunctions\inputfiles\wind.sum', 'wt');
    fprintf(fid, 'T\tCLOCKWISE\n');
    fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', H);
    fprintf(fid, '%0.3f\tUBAR\n', U);
    fprintf(fid, '%0.3f\tTI(u)\n', 100);
    fprintf(fid, '%0.3f\tTI(v)\n', 100);
    fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
    fprintf(fid, '0\tHEIGHT OFFSET');
    fclose(fid);
    
end

%% EDC (10)
if Wind.Type == 10
    
    % IEC class
    if Wind.Class(2) == 1
        Iref = 0.16;
    elseif Wind.Class(2) == 2
        Iref = 0.14;
    elseif Wind.Class(2) == 3
        Iref = 0.12;
    end
    if H >= 60
        Lambda1 = 0.7*H;
    else
        Lambda1 = 42;
    end
    
    % NTM values
    alpha = 0.2;
    sigma1 = Iref*(0.75*U+5.6);
    
    % EDC parameters
    T0 = Wind.EDC;
    Tg = 6;

    % Wind profile
    t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
    if rem(length(t),2) ~= 0
        t = [t, t+Wind.dt];
    end
    x = U*t;
    y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
    z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
    [T,~,Z] = ndgrid(t,y,z);
    thetae = 4*atan(sigma1/(U*(1+0.1*2*R/Lambda1)));
    theta = 0.5*thetae * (1-cos(pi*(T-T0)/Tg));
    theta(T < T0) = 0;
    theta(T > T0+Tg) = thetae;
    u = U*(Z/H).^alpha .* cos(theta);
    v = U*(Z/H).^alpha .* sin(theta);
    w = zeros(size(u));
    
    % Save .wnd file
    writebladed('subfunctions\inputfiles\wind',(u-U)/U,v/U,w/U,x,y,z,U);
    
    % Save .sum file
    fid = fopen('subfunctions\inputfiles\wind.sum', 'wt');
    fprintf(fid, 'T\tCLOCKWISE\n');
    fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', H);
    fprintf(fid, '%0.3f\tUBAR\n', U);
    fprintf(fid, '%0.3f\tTI(u)\n', 100);
    fprintf(fid, '%0.3f\tTI(v)\n', 100);
    fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
    fprintf(fid, '0\tHEIGHT OFFSET');
    fclose(fid);
    
end

%% ECG (11)
if Wind.Type == 11
    
    % IEC class
    if Wind.Class(2) == 1
        Iref = 0.16;
    elseif Wind.Class(2) == 2
        Iref = 0.14;
    elseif Wind.Class(2) == 3
        Iref = 0.12;
    end
    
    % NTM values
    alpha = 0.2;
    
    % ECG parameters
    T0 = Wind.ECG;
    Tg = 10;
    Vcg = 15;

    % Wind profile
    t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
    if rem(length(t),2) ~= 0
        t = [t, t+Wind.dt];
    end
    x = U*t;
    y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
    z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
    [T,~,Z] = ndgrid(t,y,z);
    V = U*(Z/H).^alpha + 0.5*Vcg*(1-cos(pi*(T-T0)/Tg));
    V(T < T0) = U*(Z(T<T0)/H).^alpha;
    V(T > T0+Tg) = U*(Z(T>T0+Tg)/H).^alpha + Vcg;
    if U < 4
        thetacg = pi;
    else
        thetacg = 4*pi/U;
    end
    theta = 0.5*thetacg*(1-cos(pi*(T-T0)/Tg));
    theta(T < T0) = 0;
    theta(T > T0+Tg) = thetacg;
    u = V .* cos(theta);
    v = V .* sin(theta);
    w = zeros(size(u));
    
    % Save .wnd file
    writebladed('subfunctions\inputfiles\wind',(u-U)/U,v/U,w/U,x,y,z,U);
    
    % Save .sum file
    fid = fopen('subfunctions\inputfiles\wind.sum', 'wt');
    fprintf(fid, 'T\tCLOCKWISE\n');
    fprintf(fid, '%0.0f\tHUB HEIGHT\n\n', H);
    fprintf(fid, '%0.3f\tUBAR\n', U);
    fprintf(fid, '%0.3f\tTI(u)\n', 100);
    fprintf(fid, '%0.3f\tTI(v)\n', 100);
    fprintf(fid, '%0.3f\tTI(w)\n\n', 100);
    fprintf(fid, '0\tHEIGHT OFFSET');
    fclose(fid);
    
end