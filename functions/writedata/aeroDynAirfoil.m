function text = aeroDynAirfoil(iFoil, Airfoil)
% aeroDynAirfoil writes data in a character array, which can be written to
% file. The writing of the character array is separated by this function
% from writing the file in FileIOClass.m, to separate agnostic fileIO from
% the generation of the content of the file.
%
%   Syntax
%     text = aeroDynAirfoil(iFoil, Airfoil)
%
%   Input arguments
%    iFoil - Index of the aerofoil in the aerofoil data set
%    Airfoil - Structure with aerofoil data
%
%   Output
%     text - Character array formatted as an AirfoilInfo v1.01.x aerofoil
%     input file for aeroDyn that can be written to a plain-text file named
%     AeroDyn_[index].dat.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%
%   Called by
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

% Zero-lift values
if max(Airfoil.Cl{iFoil}) > 0
    
    Cli = Airfoil.Cl{iFoil};
    Cli(Airfoil.Alpha{iFoil} < -45) = 1e6;
    Cli(Airfoil.Alpha{iFoil} > 45) = 1e6;
    i0 = find(abs(Cli) == min(abs(Cli)));
    i0 = i0(ceil(length(i0)/2));
    alpha_0 = Airfoil.Alpha{iFoil}(i0);
    Cd_0 = Airfoil.Cd{iFoil}(i0);
    Cm_0 = Airfoil.Cm{iFoil}(i0);
    
else
    
    alpha_0 = 0;
    Cd_0 = min(Airfoil.Cd{iFoil});
    Cm_0 = 0;
    
end

text = '';
text = append(text, sprintf('! ------------ AirfoilInfo v1.01.x Input File ----------------------------------\n'));
text = append(text, sprintf(['!', Airfoil.Name{iFoil}, ' properties\n']));
text = append(text, sprintf('!Created %s.\n', datetime));
text = append(text, sprintf('! note that this file uses Marshall Buhl''s new input file processing; start all comment lines with !\n'));
text = append(text, sprintf('! ------------------------------------------------------------------------------\n'));
text = append(text, sprintf('"DEFAULT"     InterpOrd         ! Interpolation order to use for quasi-steady table lookup {1=linear; 3=cubic spline; "default"} [default=3]\n'));
text = append(text, sprintf('          1   NonDimArea        ! The non-dimensional area of the airfoil (area/chord^2) (set to 1.0 if unsure or unneeded)\n'));
text = append(text, sprintf('          0   NumCoords         ! The number of coordinates in the airfoil shape file.  Set to zero if coordinates not included.\n'));
text = append(text, sprintf('          1   NumTabs           ! Number of airfoil tables in this file.  Each table must have lines for Re and Ctrl.\n'));
text = append(text, sprintf('! ------------------------------------------------------------------------------\n'));
text = append(text, sprintf('! data for table 1\n'));
text = append(text, sprintf('! ------------------------------------------------------------------------------\n'));
text = append(text, sprintf('       0.75   Re                ! Reynolds number in millions\n'));
text = append(text, sprintf('          0   Ctrl              ! Control setting (must be 0 for current AirfoilInfo)\n'));
text = append(text, sprintf('True          InclUAdata        ! Is unsteady aerodynamics data included in this table? If TRUE, then include 30 UA coefficients below this line\n'));
text = append(text, sprintf('!........................................\n'));
text = append(text, sprintf('      %5.4f   alpha0            ! 0-lift angle of attack, depends on airfoil.\n', alpha_0));
text = append(text, sprintf('      %5.4f   alpha1            ! Angle of attack at f=0.7, (approximately the stall angle) for AOA>alpha0. (deg)\n', Airfoil.StallAngle1(iFoil)));
text = append(text, sprintf('      %5.4f   alpha2            ! Angle of attack at f=0.7, (approximately the stall angle) for AOA<alpha0. (deg)\n', Airfoil.StallAngle2(iFoil)));
text = append(text, sprintf('          0   eta_e             ! Recovery factor in the range [0.85 - 0.95] used only for UAMOD=1, it is set to 1 in the code when flookup=True. (-)\n'));
text = append(text, sprintf('      %5.4f   C_nalpha          ! Slope of the 2D normal force coefficient curve. (1/rad)\n', Airfoil.CnSlope(iFoil)));
text = append(text, sprintf('"DEFAULT"     T_f0              ! Initial value of the time constant associated with Df in the expression of Df and f''. [default = 3]\n'));
text = append(text, sprintf('"DEFAULT"     T_V0              ! Initial value of the time constant associated with the vortex lift decay process; it is used in the expression of Cvn. It depends on Re,M, and airfoil class. [default = 6]\n'));
text = append(text, sprintf('"DEFAULT"     T_p               ! Boundary-layer,leading edge pressure gradient time constant in the expression of Dp. It should be tuned based on airfoil experimental data. [default = 1.7]\n'));
text = append(text, sprintf('"DEFAULT"     T_VL              ! Initial value of the time constant associated with the vortex advection process; it represents the non-dimensional time in semi-chords, needed for a vortex to travel from LE to trailing edge (TE)); it is used in the expression of Cvn. It depends on Re, M (weakly), and airfoil. [valid range = 6 - 13, default = 11]\n'));
text = append(text, sprintf('"DEFAULT"     b1                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.14]\n'));
text = append(text, sprintf('"DEFAULT"     b2                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.53]\n'));
text = append(text, sprintf('"DEFAULT"     b5                ! Constant in the expression of K''''_q,Cm_q^nc, and k_m,q.  [from  experimental results, defaults to 5]\n'));
text = append(text, sprintf('"DEFAULT"     A1                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.3]\n'));
text = append(text, sprintf('"DEFAULT"     A2                ! Constant in the expression of phi_alpha^c and phi_q^c.  This value is relatively insensitive for thin airfoils, but may be different for turbine airfoils. [from experimental results, defaults to 0.7]\n'));
text = append(text, sprintf('"DEFAULT"     A5                ! Constant in the expression of K''''_q,Cm_q^nc, and k_m,q. [from experimental results, defaults to 1]\n'));
text = append(text, sprintf('          0   S1                ! Constant in the f curve best-fit for alpha0<=AOA<=alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   S2                ! Constant in the f curve best-fit for         AOA> alpha1; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   S3                ! Constant in the f curve best-fit for alpha2<=AOA< alpha0; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   S4                ! Constant in the f curve best-fit for         AOA< alpha2; by definition it depends on the airfoil. [ignored if UAMod<>1]\n'));
text = append(text, sprintf('      %5.4f   Cn1               ! Critical value of C0n at leading edge separation. It should be extracted from airfoil data at a given Mach and Reynolds number. It can be calculated from the static value of Cn at either the break in the pitching moment or the loss of chord force at the onset of stall. It is close to the condition of maximum lift of the airfoil at low Mach numbers.\n', Airfoil.CritCn1(iFoil)));
text = append(text, sprintf('      %5.4f   Cn2               ! As Cn1 for negative AOAs.\n', Airfoil.CritCn2(iFoil)));
text = append(text, sprintf('       0.19   St_sh             ! Strouhal''s shedding frequency constant.  [default = 0.19]\n'));
text = append(text, sprintf('      %5.4f   Cd0               ! 2D drag coefficient value at 0-lift.\n', Cd_0));
text = append(text, sprintf('      %5.4f   Cm0               ! 2D pitching moment coefficient about 1/4-chord location, at 0-lift, positive if nose up. [If the aerodynamics coefficients table does not include a column for Cm, this needs to be set to 0.0]\n', Cm_0));
text = append(text, sprintf('          0   k0                ! Constant in the hat(x)_cp curve best-fit; = (hat(x)_AC-0.25).  [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   k1                ! Constant in the hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   k2                ! Constant in the hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   k3                ! Constant in the hat(x)_cp curve best-fit.  [ignored if UAMod<>1]\n'));
text = append(text, sprintf('          0   k1_hat            ! Constant in the expression of Cc due to leading edge vortex effects.  [ignored if UAMod<>1]\n'));
text = append(text, sprintf('"DEFAULT"     x_cp_bar          ! Constant in the expression of hat(x)_cp^v. [ignored if UAMod<>1, default = 0.2]\n'));
text = append(text, sprintf('"DEFAULT"     UACutout          ! Angle of attack above which unsteady aerodynamics are disabled (deg). [Specifying the string "Default" sets UACutout to 45 degrees]\n'));
text = append(text, sprintf('"DEFAULT"     filtCutOff        ! Cut-off frequency (-3 dB corner frequency) for low-pass filtering the AoA input to UA, as well as the 1st and 2nd derivatives (Hz) [default = 20]\n'));
text = append(text, sprintf('!........................................\n'));
text = append(text, sprintf('! Table of aerodynamics coefficients\n'));
text = append(text, sprintf('          %i   NumAlf            ! Number of data lines in the following table\n', length(Airfoil.Alpha{iFoil})));
text = append(text, sprintf('!    Alpha      Cl      Cd        Cm\n'));
text = append(text, sprintf('!    (deg)      (-)     (-)       (-)\n'));
for j = 1:length(Airfoil.Alpha{iFoil})
    text = append(text, sprintf('%5.4f    %5.4f    %5.4f    %5.4f \n', Airfoil.Alpha{iFoil}(j), Airfoil.Cl{iFoil}(j), Airfoil.Cd{iFoil}(j), Airfoil.Cm{iFoil}(j)));
end
text = append(text, sprintf('! ------------------------------------------------------------------------------\n'));
