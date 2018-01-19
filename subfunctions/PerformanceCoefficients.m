function [CT, CQ] = PerformanceCoefficients(Blade, Airfoil, pitch, lambda)

r = Blade.Radius/Blade.Radius(end);

CQ = 0;
CT = 0;

for k = 1:(length(Blade.Radius)-1)
    % Initialise search for induction factor (where blade element and
    % momentum theory give equal thrust coefficient)
    a1 = -0.4;
    a2 = 0.80; % 0.99; Induction factor can go up to 1 with this model, but a2 = 0.99 gives convergence problems
    [~, CTBE1, CTMT1] = BEM(a1, k, Blade, Airfoil, pitch, lambda);
    [~, CTBE2, CTMT2] = BEM(a2, k, Blade, Airfoil, pitch, lambda);
    Diff1 = CTBE1 - CTMT1;
    Diff2 = CTBE2 - CTMT2;
    success = false;
    
    for iter = 1:100
        a = a1 + abs(Diff1/(Diff1-Diff2))*(a2-a1);
        [dCQ, CTBE, CTMT] = BEM(a, k, Blade, Airfoil, pitch, lambda);
        Diff = CTBE - CTMT;
        if abs(Diff) < 0.005
            success = true;
            break
        end
        if Diff*Diff1 > 0 % Difference of equal sign
           a1 = a;
           Diff1 = Diff;
        else
           a2 = a;
           Diff2 = Diff;
        end
    end
    if ~success
        warning('Tolerance not met during iteration of axial induction factor')
    end
    
    % Add coefficients
    % (2*r is a weighing factor for increasing annulus with r)
    CQ = CQ+(r(k+1)-r(k))*dCQ*2*r(k);
    dCT = (r(k+1)-r(k))*(CTBE+CTMT)/2;
    CT = CT+dCT*2*r(k);
end
    
function [dCQ, CTBE, CTMT] = BEM(a, k, Blade, Airfoil, pitch, lambda)

r = Blade.Radius(k)/Blade.Radius(end);
lambdar = lambda * r;
sigmar = Blade.Number*Blade.Chord(k)/(2*pi*Blade.Radius(k));

% Iteration of tangential induction factor
a_old = 0;
success = false;
for iter = 1:100
    % Inflow angleclose all force
    phi = atan((1-a)/((1+a_old)*lambdar));

    % Tip loss correction
    F = 2/pi * acos(exp(-Blade.Number/2*(1-r)/(r*sin(abs(phi)))));

    % Angle of attack
    alpha = phi*180/pi - pitch - Blade.Twist(k);

    % Aerodynamic force coefficients
    [~,ia] = unique(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))});
    Cl = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}(ia), Airfoil.Cl{Blade.IFoil(Blade.NFoil(k))}(ia), alpha);
    Cd = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}(ia), Airfoil.Cd{Blade.IFoil(Blade.NFoil(k))}(ia), alpha);

    % Tangential induction factor
    a_ = 1/(4*F*cos(phi)/(sigmar*Cl)-1);

    % Convergence check
    if abs(a_ - a_old) < 0.002
        success = true;
        break
    end
    a_old = a_;
end
if ~success
    warning('Tolerance not met during iteration of tangential induction factor')
end

% Thrust coefficient per rotor annulus (blade element theory)
CTBE = sigmar*(1-a)^2*(Cl*cos(phi)+Cd*sin(phi))/(sin(phi)^2);
    
% Thrust coefficient per rotor annulus (momentum theory)
if a <= 0.4
    CTMT = 4*a*(1-a)*F;
else
    CTMT = 8/9 + (4*F-40/9)*a + (50/9-4*F)*a^2; % High induction (heavy loading) according to Buhl
end

% Torque coefficient per rotor annulus 
% Derived directly from tangential force
dCQ = sigmar*r*(Cl*sin(phi)-Cd*cos(phi))*(1-a)^2/(sin(phi)^2);
%{
% Derived from contribution to torque according to Hansen
% (gives almost equal results)
dCQ = lambda*(1-a)*(1+a_)*r^2*sigmar*(Cl*sin(phi)-Cd*cos(phi))/(sin(phi)*cos(phi));
%}



