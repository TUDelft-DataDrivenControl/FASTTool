function AeroDyn(Blade,Airfoil,Tower,mode,AirDensity)

if contains(mode,'Linearize')
    AFAeroMod = 1;
    WakeMod = 1;
else
    AFAeroMod = 2;
    WakeMod = 1;
end

%% AeroDyn input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'AeroDyn.dat'], 'wt');
fprintf(fid, '------- AERODYN v15.03.* INPUT FILE ------------------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '======  General Options  ============================================================================\n');
fprintf(fid, 'False         Echo               - Echo the input to "<rootname>.AD.ech"?  (flag)\n');
fprintf(fid, '"default"     DTAero             - Time interval for aerodynamic calculations {or "default"} (s)\n');
fprintf(fid, '          %i   WakeMod            - Type of wake/induction model (switch) {0=none, 1=BEMT}\n', WakeMod);
fprintf(fid, '          %i   AFAeroMod          - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model}\n', AFAeroMod);
fprintf(fid, '          0  TwrPotent          - Type tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}\n');
fprintf(fid, 'False          TwrShadow          – Calculate tower influence on wind based on downstream tower shadow? (flag)\n');
fprintf(fid, 'False           TwrAero            - Calculate tower aerodynamic loads? (flag)\n');
fprintf(fid, 'False          FrozenWake         - Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]\n');
fprintf(fid, '======  Environmental Conditions  ===================================================================\n');
fprintf(fid, '      %4.3f   AirDens            - Air density (kg/m^3)\n', AirDensity);
fprintf(fid, '  1.464E-05   KinVisc            - Kinematic air viscosity (m^2/s)\n');
fprintf(fid, '        335   SpdSound           - Speed of sound (m/s)\n');
fprintf(fid, '======  Blade-Element/Momentum Theory Options  ====================================================== [used only when WakeMod=1]\n');
fprintf(fid, '          2   SkewMod            - Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1]\n');
fprintf(fid, 'True          TipLoss            - Use the Prandtl tip-loss model? (flag) [used only when WakeMod=1]\n');
fprintf(fid, 'True          HubLoss            - Use the Prandtl hub-loss model? (flag) [used only when WakeMod=1]\n');
fprintf(fid, 'True          TanInd             - Include tangential induction in BEMT calculations? (flag) [used only when WakeMod=1]\n');
fprintf(fid, 'False         AIDrag             - Include the drag term in the axial-induction calculation? (flag) [used only when WakeMod=1]\n');
fprintf(fid, 'False         TIDrag             - Include the drag term in the tangential-induction calculation? (flag) [used only when WakeMod=1 and TanInd=TRUE]\n');
fprintf(fid, '"Default"     IndToler           - Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [used only when WakeMod=1]\n');
fprintf(fid, '        100   MaxIter            - Maximum number of iteration steps (-) [used only when WakeMod=1]\n');
fprintf(fid, '======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]\n');
fprintf(fid, '          3   UAMod              - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAeroMod=2]\n');
fprintf(fid, 'True          FLookup            - Flag to indicate whether a lookup for f will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]\n');
fprintf(fid, '======  Airfoil Information =========================================================================\n');
fprintf(fid, '          1   InCol_Alfa         - The column in the airfoil tables that contains the angle of attack (-)\n');
fprintf(fid, '          2   InCol_Cl           - The column in the airfoil tables that contains the lift coefficient (-)\n');
fprintf(fid, '          3   InCol_Cd           - The column in the airfoil tables that contains the drag coefficient (-)\n');
fprintf(fid, '          4   InCol_Cm           - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)\n');
fprintf(fid, '          0   InCol_Cpmin        - The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)\n');
fprintf(fid, '          %i   NumAFfiles         - Number of airfoil files used (-)\n', length(Blade.IFoil));
fprintf(fid, '"%s"                         AFNames            - Airfoil file names (NumAFfiles lines) (quoted strings)\n', 'AeroDyn_Cylinder 1.dat');
for i = 2:length(Blade.IFoil)
    fprintf(fid, '"%s"\n', ['AeroDyn_', Airfoil.Name{Blade.IFoil(i)}, '.dat']);
end
fprintf(fid, '======  Rotor/Blade Properties  =====================================================================\n');
fprintf(fid, 'True          UseBlCm            - Include aerodynamic pitching moment in calculations?  (flag)\n');
fprintf(fid, '"AeroDyn_blade.dat"    ADBlFile(1)        - Name of file containing distributed aerodynamic properties for Blade #1 (-)\n');
fprintf(fid, '"AeroDyn_blade.dat"    ADBlFile(2)        - Name of file containing distributed aerodynamic properties for Blade #2 (-) [unused if NumBl < 2]\n');
fprintf(fid, '"AeroDyn_blade.dat"    ADBlFile(3)        - Name of file containing distributed aerodynamic properties for Blade #3 (-) [unused if NumBl < 3]\n');
fprintf(fid, '======  Tower Influence and Aerodynamics ============================================================= [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]\n');
fprintf(fid, '          %i   NumTwrNds         - Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]\n', length(Tower.Height));
fprintf(fid, 'TwrElev        TwrDiam        TwrCd\n');
fprintf(fid, '(m)              (m)           (-)\n');
for i = 1:length(Tower.Height)
    fprintf(fid, '%5.4f    %5.4f    %5.4f\n', Tower.Height(i), Tower.Diameter(i), 0.6);
end
fprintf(fid, '======  Outputs  ====================================================================================\n');
fprintf(fid, 'False         SumPrint            - Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)\n');
fprintf(fid, '          0   NBlOuts             - Number of blade node outputs [0 - 9] (-)\n');
fprintf(fid, ' 1, 9, 19     BlOutNd             - Blade nodes whose values will be output  (-)\n');
fprintf(fid, '          0   NTwOuts             - Number of tower node outputs [0 - 9]  (-)\n');
fprintf(fid, ' 1, 2, 3, 4, 5     TwOutNd             - Tower nodes whose values will be output  (-)\n');
fprintf(fid, '                   OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)\n');
fprintf(fid, 'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n');
fprintf(fid, '---------------------------------------------------------------------------------------\n');
fclose(fid);

%% Blade input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'AeroDyn_blade.dat'], 'wt');
fprintf(fid, '------- AERODYN v15.00.* BLADE DEFINITION INPUT FILE -------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '======  Blade Properties =================================================================\n');
fprintf(fid, '         %i   NumBlNds           - Number of blade nodes used in the analysis (-)\n', length(Blade.Radius));
fprintf(fid, '  BlSpn        BlCrvAC        BlSwpAC        BlCrvAng       BlTwist        BlChord          BlAFID\n');
fprintf(fid, '   (m)           (m)            (m)            (deg)         (deg)           (m)              (-)\n');
Blade.Radius(end) = Blade.Radius(end) - 1e-4;
for i = 1:length(Blade.Radius)
    fprintf(fid, '%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\t%7.7E\t%i\n', ...
        Blade.Radius(i)-Blade.Radius(1), ...
        0, ...
        0, ...
        0, ...
        Blade.Twist(i), ...
        Blade.Chord(i), ...
        Blade.NFoil(i));
end
fclose(fid);

%% Airfoil input files
for i = 1:length(Blade.IFoil)

    % Zero-lift values
    if max(Airfoil.Cl{Blade.IFoil(i)}) > 0
        
        Cli = Airfoil.Cl{Blade.IFoil(i)};
        Cli(Airfoil.Alpha{Blade.IFoil(i)} < -45) = 1e6;
        Cli(Airfoil.Alpha{Blade.IFoil(i)} > 45) = 1e6;
        i0 = find(abs(Cli) == min(abs(Cli)));
        i0 = i0(ceil(length(i0)/2));
        alpha_0 = Airfoil.Alpha{Blade.IFoil(i)}(i0);
        Cd_0 = Airfoil.Cd{Blade.IFoil(i)}(i0);
        Cm_0 = Airfoil.Cm{Blade.IFoil(i)}(i0);
        
    else
        
        alpha_0 = 0;
        Cd_0 = min(Airfoil.Cd{Blade.IFoil(i)});
        Cm_0 = 0;
        
    end

    fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'AeroDyn_', Airfoil.Name{Blade.IFoil(i)}, '.dat'], 'wt');
    fprintf(fid, '! ------------ AirfoilInfo v1.01.x Input File ----------------------------------\n');
    fprintf(fid, ['!', Airfoil.Name{Blade.IFoil(i)}, ' properties\n']);
    fprintf(fid, '!Created %s.\n', datestr(now));
    fprintf(fid, '! note that this file uses Marshall Buhl''s new input file processing; start all comment lines with !\n');
    fprintf(fid, '! ------------------------------------------------------------------------------\n');
    fprintf(fid, '"DEFAULT"     InterpOrd         ! Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]\n');
    fprintf(fid, '          1   NonDimArea        ! The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)\n');
    fprintf(fid, '          0   NumCoords         ! The number of coordinates in the airfoil shape file.  Set to zero if coordinates not included.\n');
    fprintf(fid, '          1   NumTabs           ! Number of airfoil tables in this file.  Each table must have lines for Re and Ctrl.\n');
    fprintf(fid, '! ------------------------------------------------------------------------------\n');
    fprintf(fid, '! data for table 1\n');
    fprintf(fid, '! ------------------------------------------------------------------------------\n');
    fprintf(fid, '       0.75   Re                ! Reynolds number in millions\n');
    fprintf(fid, '          0   Ctrl              ! Control setting (must be 0 for current AirfoilInfo)\n');
    fprintf(fid, 'True          InclUAdata        ! Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line\n');
    fprintf(fid, '!........................................\n');
    fprintf(fid, '      %5.4f   alpha0            ! 0-lift angle of attack, depends on airfoil.\n', alpha_0);
    fprintf(fid, '      %5.4f   alpha1            ! Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)\n', Airfoil.StallAngle1(Blade.IFoil(i)));
    fprintf(fid, '      %5.4f   alpha2            ! Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)\n', Airfoil.StallAngle2(Blade.IFoil(i)));
    fprintf(fid, '          0   eta_e             ! Recovery factor in the range [0.85 - 0.95] used only for UAMOD=1, it is set to 1 in the code when flookup=True. (-)\n');
    fprintf(fid, '      %5.4f   C_nalpha          ! Slope of the 2D normal force coefficient curve. (1/rad)\n', Airfoil.CnSlope(Blade.IFoil(i)));
    fprintf(fid, '"DEFAULT"     T_f0              ! Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]\n');
    fprintf(fid, '"DEFAULT"     T_V0              ! Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]\n');
    fprintf(fid, '"DEFAULT"     T_p               ! Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]\n');
    fprintf(fid, '"DEFAULT"     T_VL              ! Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]\n');
    fprintf(fid, '"DEFAULT"     b1                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]\n');
    fprintf(fid, '"DEFAULT"     b2                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]\n');
    fprintf(fid, '"DEFAULT"     b5                ! Constant in the expression of K''''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]\n');
    fprintf(fid, '"DEFAULT"     A1                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]\n');
    fprintf(fid, '"DEFAULT"     A2                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]\n');
    fprintf(fid, '"DEFAULT"     A5                ! Constant in the expression of K''''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]\n');
    fprintf(fid, '          0   S1                ! Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   S2                ! Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   S3                ! Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   S4                ! Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UAMod<>1]\n');
    fprintf(fid, '      %5.4f   Cn1               ! Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.\n', Airfoil.CritCn1(Blade.IFoil(i)));
    fprintf(fid, '      %5.4f   Cn2               ! As Cn1 for negative AOAs.\n', Airfoil.CritCn2(Blade.IFoil(i)));
    fprintf(fid, '       0.19   St_sh             ! Strouhal''s shedding frequency constant.  [default = 0.19]\n');
    fprintf(fid, '      %5.4f   Cd0               ! 2D drag coefficient value at 0-lift.\n', Cd_0);
    fprintf(fid, '      %5.4f   Cm0               ! 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]\n', Cm_0);
    fprintf(fid, '          0   k0                ! Constant in the hat(x)_cp curve best-fit; = (hat(x)_AC-0.25).  [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   k1                ! Constant in the hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   k2                ! Constant in the hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   k3                ! Constant in the hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n');
    fprintf(fid, '          0   k1_hat            ! Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UAMod<>1]\n');
    fprintf(fid, '"DEFAULT"     x_cp_bar          ! Constant in the expression of hat(x)_cp^v. [ignored if UAMod<>1, default = 0.2]\n');
    fprintf(fid, '"DEFAULT"     UACutout          ! Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UACutout to 45 degrees]\n');
    fprintf(fid, '"DEFAULT"     filtCutOff        ! Cut-off frequency (-3 dB corner frequency) for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (Hz) [default = 20]\n');
    fprintf(fid, '!........................................\n');
    fprintf(fid, '! Table of aerodynamics coefficients\n');
    fprintf(fid, '          %i   NumAlf            ! Number of data lines in the following table\n', length(Airfoil.Alpha{Blade.IFoil(i)}));
    fprintf(fid, '!    Alpha      Cl      Cd        Cm\n');
    fprintf(fid, '!    (deg)      (-)     (-)       (-)\n');
    for j = 1:length(Airfoil.Alpha{Blade.IFoil(i)})
        fprintf(fid, '%5.4f    %5.4f    %5.4f    %5.4f \n', Airfoil.Alpha{Blade.IFoil(i)}(j), Airfoil.Cl{Blade.IFoil(i)}(j), Airfoil.Cd{Blade.IFoil(i)}(j), Airfoil.Cm{Blade.IFoil(i)}(j));
    end
    fprintf(fid, '! ------------------------------------------------------------------------------\n');
    fclose(fid);
end
