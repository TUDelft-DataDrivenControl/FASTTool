function [CT, CQ, OmegaU, PitchAngle] = steadyState(Blade, Airfoil, Drivetrain, Control, wind, AirDensity)
% steadyState determines the thrust and torque coefficients, the rotational
% speed and the pitch angle as a function of wind speed
%
%   Syntax
%     [CT, CQ, OmegaU, PitchAngle] = steadyState(Blade, Airfoil, Drivetrain, Control, wind, AirDensity)
%
%   Input arguments
%     Blade - Structure with blade design information
%     Airfoil - Structure with aerofoil data
%     pitch - Pitch angle at which performance needs to be determined
%     lambda - Tip speed ratio at which performance needs to be determined
%
%   Output
%     CT - Thrust coefficients for all wind speeds
%     CQ - Torque coefficients for all wind speeds
%     OmegaU - Rotational speeds of the rotor for all wind speeds
%     PitchAngle - Pitch angles for all wind speeds
%
%   Called by
%     WindTurbineClass.performSteadyOpAnalysis
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

% The operational conditions for each wind speed are determined by looking
% for the point where the rotor-performance curve of torque versus
% rotational speed crosses the torque-control curve (which is expressed in
% the same two parameters). Calculations are done for torque and rotational
% speed in the high-speed shaft.

% Tracing information of warnings is suppressed, to avoid that too much
% information is output to the command window.
warnStruct = warning('off', 'backtrace');

% Initialisation of each parameter or group of parameters will be explained
% separately.

% Output arrays are initialysed to give them the right dimensions for use
% in the loop over wind speeds. The initialisation values (zeros) remain
% for wind speeds that fall outside the operational range from cut-in to
% cut-out.
OmegaU = zeros(size(wind));
PitchAngle = Control.Pitch.Max*ones(size(wind));
CT = zeros(size(wind));
CQ = zeros(size(wind));

% Torque control speeds are copied and convert to rad/s for convenience of
% calculations
OmegaA = Control.Torque.SpeedA*2*pi/60;
OmegaB = Control.Torque.SpeedB*2*pi/60;
OmegaB2 = Control.Torque.SpeedB2*2*pi/60;
OmegaC = Control.Torque.SpeedC*2*pi/60;

% The indication of the control region is initialysed. The control region
% can be set to:
%
% 'lin' = linear control curve
% 'square' = quadratic control curve
% 'point' = single control point
%
% For different control regions, different algorithms are used to determine
% the operational conditions.
%
% It is initialised with 'lin', because when rotor torque is initialised
% with 0, the analysis first encounters the linear control region from A
% (minimum operations speed) to B (rotational speed at which the optimal
% mode gain is obtained).
Region = 'lin';

% The crossing of the rotor-performance and torque-control curves is done
% for segments that are defined by two points of the respective curves,
% over the same range for rotational speeds. The rotational speed range is
% initialised to the first part of the control curve, which is defined by
% [OmegaA, OmegaB]. The corresponding torques for the control curve are
% determined from the control settings. The rotor-performance is set to
% zero over the entire range, since it will be updated as first thing in
% the loop over (operational) wind speeds. The blade pitch angles at the
% bounds of the range of rotational speeds are initialised at the
% fine-pitch angle, since these speeds are in the partial load region.
% There may be no crossing between the two initialysed curves (with the
% updated torque for the rotor performance). How the algorithm updates the
% curves to find a crossing is explained at the start of the loop over wind
% speeds, after the update of the rotor torque curve. 
Omega = [OmegaA, OmegaB]; % Rotational speeds
Qc = [0, Control.Torque.OptGain*OmegaB^2]; % Torque-control curve
Qr = [0, 0]; % Rotor performance curve
Beta = [Control.Pitch.Fine, Control.Pitch.Fine]; % Pitch angle

% An boolean variable is used to identify the transition from the
% partial-load region ('lin') to the full-load region ('square'). The first
% wind speed in the region 'square' requires a different approach to
% search for the crossing of the curves than later points.
FirstOptimal = true;

% A tolerance for convergence is set on the difference in torque from the
% control curve and from the rotor-performance curve, relative to demanded
% torque.
Tolerance = 0.0005;

for i = 1:length(wind)
    % Initialysed values of the output parameters are only changed for wind
    % speeds inside the operational range.
    if wind(i) >= Control.WindSpeed.Cutin && wind(i) <= Control.WindSpeed.Cutout

        % Update of rotor torque range for the current wind speed for the
        % current rotational speed range and pitch angles.
        [~, CQr] = performanceCoefficients(Blade, Airfoil, Beta(1), (Omega(1)/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/wind(i));
        Qr(1) = 0.5*CQr*AirDensity*wind(i)^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
        [~, CQr] = performanceCoefficients(Blade, Airfoil, Beta(2), (Omega(2)/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/wind(i));
        Qr(2) = 0.5*CQr*AirDensity*wind(i)^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;

        % The next code checks whether for the current rotational speed
        % range there is a crossing between the torque-control curve and
        % the rotor-torque curve. If not, it increases the rotational speed
        % range to the next control region and tests again. This continues,
        % until the rated rotational speed, OmegaC, is reached. The outcome
        % of this process is the identification of the control region and
        % corresponding rotational speed range in which to search for the
        % operational point. For control regions 'lin' and 'square', the
        % procedure ensures that a rotational speed can be found for which
        % torque on the control curve equals torque on the
        % rotor-performance curve. For control region 'point', which marks
        % full load operation, the algorithm sets the rated speed, for
        % which the pitch angle can be found at which torque on the control
        % curve equals torque on the rotor-performance curve.
        % These searches for corresponding torque on both curves are
        % performed later, in the code starting with 'switch Region'.
        while (Omega(1) < OmegaC && Qc(2) < Qr(2))

            Omega(1) = Omega(2);
            Qc(1) = Qc(2);
            Qr(1) = Qr(2);

            % Determine control region
            if Omega(1) == OmegaB % Region 2
                Omega(2) = OmegaB2;
                Qc(2) = Control.Torque.OptGain*OmegaB2^2;
                Region = 'square';
            elseif Omega(1) == OmegaB2 % Region 2 1/2
                Omega(2) = OmegaC;
                Qc(2) = Control.Torque.Demanded;
                Region = 'lin';
            else % Region 3
                Region = 'point';
            end

            % Determine torque from rotor performance
            [~, CQr] = performanceCoefficients(Blade, Airfoil, Beta(2), (Omega(2)/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/wind(i));
            Qr(2) = 0.5*CQr*AirDensity*wind(i)^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
        end

        % The next code searched for the operational point, at which torque
        % from the rotor-performance curve equals torque on the control
        % curve. In the cases 'lin' and 'square', this is done by searching
        % for the rotational speed at which this occurs. In this search
        % the pitch angle is left unchanged, at the fine-pitch angle. In
        % the case 'point', this is done by searching for the pitch angle
        % at which this occurs. In this search the rotational speed is left
        % unchanged, at the rated speed.
        switch Region
            case 'lin'
                Omega1 = Omega(1);
                Omega2 = Omega(2);
                Qr1 = Qr(1);
                Qr2 = Qr(2);
                x1 = Omega(1);
                y1 = Qc(1);
                x2 = Omega(2);
                y2 = Qc(2);
                success = false;

                for iter = 1:20
                    % Find crossing point of control curve and rotor torque curve
                    % http://www.ambrsoft.com/MathCalc/Line/TwoLinesIntersection/TwoLinesIntersection.htm
                    x3 = Omega1;
                    y3 = Qr1;
                    x4 = Omega2;
                    y4 = Qr2;

                    OmegaCross = ((x2*y1-x1*y2)*(x4-x3)-(x4*y3-x3*y4)*(x2-x1))/...
                        ((x2-x1)*(y4-y3)-(x4-x3)*(y2-y1));
                    QcCross = Qc(1) + (OmegaCross-Omega(1))*(Qc(2)-Qc(1))/(Omega(2)-Omega(1));
                    % Determine torque from rotor performance
                    [CTr, CQr] = performanceCoefficients(Blade, Airfoil, Beta(2), (OmegaCross/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/wind(i));
                    QrCross = 0.5*CQr*AirDensity*wind(i)^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
                    Diff = (QrCross - QcCross)/Control.Torque.Demanded;

                    if abs(Diff) < Tolerance
                        success = true;
                        break;
                    end

                    if Diff > 0 % Rotor torque at iterated crossing point lies above control curve
                        Omega1 = OmegaCross;
                        Qr1 = QrCross;
                    else
                        Omega2 = OmegaCross;
                        Qr2 = QrCross;
                    end
                end
                if ~success
                    warning('Tolerance not met during iteration of rotor speed');
                end

                OmegaU(i) = OmegaCross/Drivetrain.Gearbox.Ratio;
                PitchAngle(i) = Control.Pitch.Fine;
                CQ(i) = CQr;
                CT(i) = CTr;

            case 'square'
                % The first operational point in this control region
                % requires a search for the crossing of the control curve
                % and the rotor performance curve. For all subsequent
                % operational points in this control region the conditions
                % can be determined from this first point, since the region
                % is characterised by constant (fine) pitch and constant
                % tip speed ratio. Therefore, also the thrust and torque
                % coefficients remain constant.
                if FirstOptimal

                    Omega1 = Omega(1);
                    Omega2 = Omega(2);
                    Qr1 = Qr(1);
                    Qr2 = Qr(2);
                    success = false;

                    for iter = 1:20
                        % Find crossing point of control curve and rotor torque
                        % curve with abc-formula
                        aa = Control.Torque.OptGain;
                        bb = -(Qr2-Qr1)/(Omega2-Omega1);
                        cc = -((Omega2*Qr1-Omega1*Qr2)/(Omega2-Omega1));

                        OmegaCross = (-bb + sqrt(bb^2-4*aa*cc))/(2*aa);
                        QcCross = Control.Torque.OptGain*OmegaCross^2;

                        % Determine torque from rotor performance
                        [CTr, CQr] = performanceCoefficients(Blade, Airfoil, Beta(2), (OmegaCross/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/wind(i));
                        QrCross = 0.5*CQr*AirDensity*wind(i)^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
                        Diff = (QrCross - QcCross)/Control.Torque.Demanded;

                        if abs(Diff) < Tolerance
                            success = true;
                            break;
                        end

                        if Diff > 0 % Rotor torque at iterated crossing point lies above control curve
                            Omega1 = OmegaCross;
                            Qr1 = QrCross;
                        else
                            Omega2 = OmegaCross;
                            Qr2 = QrCross;
                        end
                    end
                    if ~success
                        warning('Tolerance not met during iteration of rotor speed');
                    end

                    OmegaU(i) = OmegaCross/Drivetrain.Gearbox.Ratio;
                    PitchAngle(i) = Control.Pitch.Fine;
                    CQ(i) = CQr;
                    CT(i) = CTr;

                    FirstOptimal = false;
                else
                    PitchAngle(i) = Control.Pitch.Fine;
                    OmegaU(i) = OmegaU(i-1) * wind(i) / wind(i-1);
                    CQ(i) = CQ(i-1);
                    CT(i) = CT(i-1);
                end

            case 'point'

                TSRFull = (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/wind(i);
                Omega(1) = OmegaC;

                % Pitch blade until the rotor torque drops below the demanded
                % torque
                while Qr(2) > Control.Torque.Demanded

                    Qr(1) = Qr(2);
                    Beta(1) = Beta(2);
                    Beta(2) = Beta(1)+0.5;

                    % Determine torque from rotor performance
                    [~, CQr] = performanceCoefficients(Blade, Airfoil, Beta(2), TSRFull);
                    Qr(2) = 0.5*CQr*AirDensity*wind(i)^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;

                end

                OmegaU(i) = OmegaC/Drivetrain.Gearbox.Ratio;
                PitchAngle(i) = interp1(Qr, Beta, Control.Torque.Demanded);
                [CTr, CQr] = performanceCoefficients(Blade, Airfoil, PitchAngle(i), TSRFull);
                CQ(i) = CQr;
                CT(i) = CTr;
        end
    end
end

warning(warnStruct); % Restore backtracing of warnings to original state

