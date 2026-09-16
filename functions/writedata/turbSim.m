function text = turbSim(Wind, U, H, windSeed)
% turbSim writes data in a character array, which can be written to file.
% The writing of the character array is separated by this function from
% writing the file in FileIOClass.m, to separate agnostic fileIO from the
% generation of the content of the file.
%
%   Syntax
%     text = turbSim(Wind, U, H, windSeed)
%
%   Input arguments
%     Wind - Structure containing specifications of the wind conditions
%     U - 10-Minute average wind speed
%     H - Hub height
%     windSeed - First seed number for the stochastic number generator
%
%   Output
%     text - Character array formatted as a TurbSim v1.06.00 input file
%     that can be written to a plain-text file named wind.inp.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%
%   Called by
%     WindTurbineClass.simulate

% TurbSim for NTM (4), EWM1 (5), EWM50 (6), ETM (8)

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

% Write TurbSim input file
text = '';
text = append(text, sprintf('TurbSim Input File. Valid for TurbSim v1.06.00, 21-Sep-2012\n'));
text = append(text, newline);
text = append(text, sprintf('---------Runtime Options-----------------------------------\n'));
text = append(text, sprintf('%i                  RandSeed1       - First random seed  (-2147483648 to 2147483647) \n', windSeed));
text = append(text, sprintf('RANLUX              RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"\n'));
text = append(text, sprintf('False               WrBHHTP         - Output hub-height turbulence parameters in binary form?  (Generates RootName.bin)\n'));
text = append(text, sprintf('False               WrFHHTP         - Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat)\n'));
text = append(text, sprintf('False               WrADHH          - Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh)\n'));
text = append(text, sprintf('False               WrADFF          - Output full-field time-series data in TurbSim/AeroDyn form? (Generates Rootname.bts)\n'));
text = append(text, sprintf('True                WrBLFF          - Output full-field time-series data in BLADED/AeroDyn form?  (Generates RootName.wnd)\n'));
text = append(text, sprintf('False               WrADTWR         - Output tower time-series data? (Generates RootName.twr)\n'));
text = append(text, sprintf('False               WrFMTFF         - Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, RootName.v, RootName.w)\n'));
text = append(text, sprintf('True                WrACT           - Output coherent turbulence time steps in AeroDyn form? (Generates RootName.cts)\n'));
text = append(text, sprintf('True                Clockwise       - Clockwise rotation looking downwind? (used only for full-field binary files - not necessary for AeroDyn)\n'));
text = append(text, sprintf('0                   ScaleIEC        - Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale uniformly; 2=use individual scales]\n'));
text = append(text, sprintf(' \n'));
text = append(text, sprintf('--------Turbine/Model Specifications-----------------------\n'));
text = append(text, sprintf('%i                  NumGrid_Z       - Vertical grid-point matrix dimension\n', Wind.Nz));
text = append(text, sprintf('%i                  NumGrid_Y       - Horizontal grid-point matrix dimension\n', Wind.Ny));
text = append(text, sprintf('%2.5f               TimeStep        - Time step [seconds]\n', Wind.dt));
text = append(text, sprintf('%2.3f               AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )\n', Wind.T));
text = append(text, sprintf('%2.3f               UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds)\n', Wind.T));
text = append(text, sprintf('%2.3f               HubHt           - Hub height [m] (should be > 0.5*GridHeight)\n', H));
text = append(text, sprintf('%2.3f               GridHeight      - Grid height [m] \n', Wind.Lz));
text = append(text, sprintf('%2.3f               GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))\n', Wind.Ly));
text = append(text, sprintf('0                   VFlowAng        - Vertical mean flow (uptilt) angle [degrees]\n'));
text = append(text, sprintf('0                   HFlowAng        - Horizontal mean flow (skew) angle [degrees]\n'));
text = append(text, sprintf('  \n'));
text = append(text, sprintf('--------Meteorological Boundary Conditions-------------------\n'));
text = append(text, sprintf('"IECKAI"            TurbModel       - Turbulence model ("IECKAI"=Kaimal, "IECVKM"=von Karman, "GP_LLJ", "NWTCUP", "SMOOTH", "WF_UPW", "WF_07D", "WF_14D", "TIDAL", or "NONE")\n'));
text = append(text, sprintf('"1-ED3"             IECstandard     - Number of IEC 61400-x standard (x=1,2, or 3 with optional 61400-1 edition number (i.e. "1-Ed2") )\n'));
text = append(text, sprintf('"%s"                IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP model, not used for other models)\n', IECturbc));
text = append(text, sprintf('"%s"                IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)\n', IEC_WindType));
text = append(text, sprintf('default             ETMc            - IEC Extreme Turbulence Model "c" parameter [m/s]\n'));
text = append(text, sprintf('"PL"                WindProfileType - Wind profile type ("JET";"LOG"=logarithmic;"PL"=power law;"H2L"=Log law for TIDAL spectral model;"IEC"=PL on rotor disk, LOG elsewhere; or "default")\n'));
text = append(text, sprintf('%2.3f               RefHt           - Height of the reference wind speed [m]\n', H));
text = append(text, sprintf('%2.3f               URef            - Mean (total) wind speed at the reference height [m/s] (or "default" for JET wind profile)\n', U));
text = append(text, sprintf('default             ZJetMax         - Jet height [m] (used only for JET wind profile, valid 70-490 m)\n'));
text = append(text, sprintf('default             PLExp           - Power law exponent [-] (or "default")           \n'));
text = append(text, sprintf('default             Z0              - Surface roughness length [m] (or "default")\n'));
text = append(text, newline);
text = append(text, sprintf('--------Non-IEC Meteorological Boundary Conditions------------\n'));
text = append(text, sprintf('default             Latitude        - Site latitude [degrees] (or "default")\n'));
text = append(text, sprintf('0.05                RICH_NO         - Gradient Richardson number \n'));
text = append(text, sprintf('default             UStar           - Friction or shear velocity [m/s] (or "default")\n'));
text = append(text, sprintf('default             ZI              - Mixing layer depth [m] (or "default")\n'));
text = append(text, sprintf('default             PC_UW           - Hub mean u''w'' Reynolds stress (or "default")\n'));
text = append(text, sprintf('default             PC_UV           - Hub mean u''v'' Reynolds stress (or "default")\n'));
text = append(text, sprintf('default             PC_VW           - Hub mean v''w'' Reynolds stress (or "default")\n'));
text = append(text, sprintf('default             IncDec1         - u-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'));
text = append(text, sprintf('default             IncDec2         - v-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'));
text = append(text, sprintf('default             IncDec3         - w-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")\n'));
text = append(text, sprintf('default             CohExp          - Coherence exponent (or "default")\n'));
text = append(text, newline);
text = append(text, sprintf('--------Coherent Turbulence Scaling Parameters-------------------\n'));
text = append(text, sprintf('"dummy"             CTEventPath     - Name of the path where event data files are located\n'));
text = append(text, sprintf('"Random"            CTEventFile     - Type of event files ("LES", "DNS", or "RANDOM")\n'));
text = append(text, sprintf('true                Randomize       - Randomize the disturbance scale and locations? (true/false)\n'));
text = append(text, sprintf('1.0                 DistScl         - Disturbance scale (ratio of wave height to rotor disk). (Ignored when Randomize = true.)\n'));
text = append(text, sprintf('0.5                 CTLy            - Fractional location of tower centerline from right (looking downwind) to left side of the dataset. (Ignored when Randomize = true.)\n'));
text = append(text, sprintf('0.5                 CTLz            - Fractional location of hub height from the bottom of the dataset. (Ignored when Randomize = true.)\n'));
text = append(text, sprintf('30.0                CTStartTime     - Minimum start time for coherent structures in RootName.cts [seconds]\n'));
text = append(text, newline);
text = append(text, sprintf('==================================================\n'));
text = append(text, sprintf('NOTE: Do not add or remove any lines in this file!\n'));
text = append(text, sprintf('==================================================\n'));
    
