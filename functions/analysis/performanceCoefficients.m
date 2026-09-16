function [CT, CQ] = performanceCoefficients(Blade, Airfoil, pitch, lambda)
% performanceCoefficients determines the thrust and torque coefficients of
% the rotor that is specified by the design in Blade and Airfoil, when
% operated at pitch angle pitch and tip speed ratio lambda
%
%   Syntax
%     [CT, CQ] = performanceCoefficients(Blade, Airfoil, pitch, lambda)
%
%   Input arguments
%     Blade - Structure with blade design information
%     Airfoil - Structure with aerofoil data
%     pitch - Pitch angle at which performance needs to be determined
%     lambda - Tip speed ratio at which performance needs to be determined
%
%   Output
%     CT - Thrust coefficient
%     CQ - Torque coefficient
%
%   Called by
%     WindTurbineClass.performRotorAnalysis
%     WindTurbineClass.linearise
%     SteadyState
%
%   Note
%     Power coefficient cP is not determined, since it can be computed from
%     cP = CQ * lambda

% Sources and equations are given in document:
% BEM - in performanceCoefficients.m.docx

% Tracing information of warnings is suppressed, to avoid that too much
% information is output to the command window.
warnStruct = warning('off', 'backtrace');

% Calculations are done as a function of dimensionless radial position r.
r = Blade.Radius/Blade.Radius(end);

% Torque and thrust coefficients of the entire blade are initialised.
% Contributions of blade elements are added to this in the for loop of k
% over the blade elements.
CQ = 0;
CT = 0;
for k = 1:(length(Blade.Radius)-1)
    % For each radial position, the induction factor is searched through
    % an iterative process. The induction factor is obtained when blade
    % element and momentum theory give an equal thrust coefficient.
    %
    % Since the difference between the thrust coefficients according to
    % blade element theory and momentum theory, respectively, is known to
    % be a monotonous function of the induction factor, and since it is
    % known to change sign over the range of the minimum and maximum
    % possible values for the induction factor, regula falsi is used as
    % root-finding algorithm for the iteration.
    %
    % The lowest value of the search range is set at a negative value,
    % since part of the blade can be in 'propeller' mode. The induction
    % factor could go up to 1 with this model, but a2 = 0.99 gave
    % convergence problems. Therefore, the highest value of the range is
    % set at 0.80.
    a1 = -0.4;
    a2 = 0.80;
    [~, CTBE1, CTMT1, ~] = BEM(a1, k, Blade, Airfoil, pitch, lambda);
    [~, CTBE2, CTMT2, ~] = BEM(a2, k, Blade, Airfoil, pitch, lambda);
    Diff1 = CTBE1 - CTMT1;
    Diff2 = CTBE2 - CTMT2;
    success = false;
    
    for iter = 1:100
        success2 = true;
        a = a1 + abs(Diff1/(Diff1-Diff2))*(a2-a1);
        [dCQ, CTBE, CTMT, BEMwarning] = BEM(a, k, Blade, Airfoil, pitch, lambda);
        if BEMwarning
            success2 = false;
        end
        Diff = CTBE - CTMT;
        if abs(Diff) < 0.005
            success = true;
            break;
        end
        if Diff*Diff1 > 0 % This means that the difference is of equal sign
           a1 = a;
           Diff1 = Diff;
        else
           a2 = a;
           Diff2 = Diff;
        end
    end
    if ~success2
        warning(['Tolerance not met during iteration of tangential induction factor for radial position ' num2str(Blade.Radius(k)) ', pitch angle ', num2str(pitch), ' deg and tip speed ratio ', num2str(lambda), '.']);
    end
    if ~success
        warning(['Tolerance not met during iteration of axial induction factor for radial position ' num2str(Blade.Radius(k)) ', pitch angle ', num2str(pitch), ' deg and tip speed ratio ', num2str(lambda), '.']);
    end
    
    % Coefficients obtained for the blade element are added to the
    % coefficients of the entire blade, as processed up to index k.
    % 2 * r is a weighing factor for increasing annulus with r.
    CQ = CQ+(r(k+1)-r(k))*dCQ*2*r(k);
    dCT = (r(k+1)-r(k))*(CTBE+CTMT)/2;
    CT = CT+dCT*2*r(k);
end

warning(warnStruct); % Restore backtracing of warnings to original state

    
function [dCQ, CTBE, CTMT, BEMwarning] = BEM(a, k, Blade, Airfoil, pitch, lambda)
% BEM determines the thrust and torque coefficients of a blade element
% according to blade element theory and its thrust coefficient according to
% momentum theory.
%
%   Syntax
%     [dCQ, CTBE, CTMT, BEMwarning] = BEM(a, k, Blade, Airfoil, pitch, lambda)
%
%   Input arguments
%     a - Induction factor
%     k - Index of blade element in Blade data
%     Blade - Structure with blade design information
%     Airfoil - Structure with aerofoil data
%     pitch - Pitch angle at which performance needs to be determined
%     lambda - Tip speed ratio at which performance needs to be determined
%
%   Output
%     dCQ - Torque coefficient of blade element
%     CTBE - Thrust coefficient of blade element according to blade element
%     theory
%     CTMT - Thrust coefficient of blade element according to momentum
%     theory
%     BEMwarning - Flag to indicate that tolerance was not met during
%     iteration of tangential induction factor
%
%   Called by
%     performanceCoefficients

% Calculations are done as a function of dimensionless radial position r.
r = Blade.Radius(k)/Blade.Radius(end);

lambdar = lambda * r;
sigmar = Blade.Number*Blade.Chord(k)/(2*pi*Blade.Radius(k));

% The tangential induction factor is searched through an iterative process.
% The inflow angle and angle of attack are determined for the latest value
% of the tangential induction factor. Subsequently, the tangential
% induction factor is updated for these angles and corresponding lift
% coefficient. This process is repeated, until the tangential induction
% factor has converged.
a_old = 0;
BEMwarning = true;
[~,ia] = unique(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}); % Subset of indices to ensure that no repetitions of alpha in the aerofoil data occur in the interpolation with interp1
for iter = 1:100
    % Determine inflow angle and angle of attack for the old tangential
    % induction factor.
    phi = atan((1-a)/((1+a_old)*lambdar)); % Inflow angle
    alpha = phi*180/pi - pitch - Blade.Twist(k); % Angle of attack 

    % Update the tangential induction factor
    Cl = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}(ia), Airfoil.Cl{Blade.IFoil(Blade.NFoil(k))}(ia), alpha); % Lift coefficient
    F = 2/pi * acos(exp(-Blade.Number/2*(1-r)/(r*sin(abs(phi))))); % Tip loss correction
    a_ = 1/(4*F*cos(phi)/(sigmar*Cl)-1);

    % Perform convergence check
    if abs(a_ - a_old) < 0.002
        Cd = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}(ia), Airfoil.Cd{Blade.IFoil(Blade.NFoil(k))}(ia), alpha); % Drag coefficient
        BEMwarning = false;
        break;
    end
    a_old = a_;
end
if BEMwarning
    Cd = max(Airfoil.Cd{Blade.IFoil(Blade.NFoil(k))}(ia)); % Drag coefficient defaults to maximum in the data for this aerofoil
end

% Determine the thrust coefficient for the annulus (blade element theory)
CTBE = sigmar*(1-a)^2*(Cl*cos(phi)+Cd*sin(phi))/(sin(phi)^2);
    
% Determine the thrust coefficient for the annulus (momentum theory)
if a <= 0.4
    % Momentum theory with tip loss correction
    CTMT = 4*a*(1-a)*F;
else
    % Momentum theory for high induction (heavy loading) according to Buhl
    CTMT = 8/9 + (4*F-40/9)*a + (50/9-4*F)*a^2;
end

% Determine the torque coefficient for the annulus. This is derived
% directly from the tangential force on the blade element (and is therefore
% based on blade element theory)
% Alternatively, it could be derived from contribution to torque according
% to Hansen:
% dCQ = lambda*(1-a)*(1+a_)*r^2*sigmar*(Cl*sin(phi)-Cd*cos(phi))/(sin(phi)*cos(phi));
% That gives almost equal results.
dCQ = sigmar*r*(Cl*sin(phi)-Cd*cos(phi))*(1-a)^2/(sin(phi)^2);


