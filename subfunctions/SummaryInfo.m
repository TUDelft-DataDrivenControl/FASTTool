%% Initialization code
function varargout = SummaryInfo(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SummaryInfo_OpeningFcn, ...
                   'gui_OutputFcn',  @SummaryInfo_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

%% Opening function
function SummaryInfo_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
Blade = varargin{1};
Airfoil = varargin{2};
Tower = varargin{3};
Nacelle = varargin{4};
Drivetrain = varargin{5};
Control = varargin{6};

% Find rated wind speed
disp('Gathering data...')
TSR = 0:0.1:20;
CP = zeros(size(TSR));
for i = 2:length(TSR)

    % Local tip speed ratio and solidity
    r = Blade.Radius/Blade.Radius(end);
    lambdar = TSR(i) * r;
    sigmar = Blade.Number*Blade.Chord./(2*pi*Blade.Radius);

    % Initial induction factors
    anew = 1/3*ones(size(Blade.Radius));
    a_new = zeros(size(Blade.Radius));
    for iter = 1:100

        a = anew;
        a_ = a_new;

        % Inflow angle
        phi = real(atan((1-a)./((1+a_).*lambdar)));

        % Tip loss correction
        F = 2/pi * acos(exp(-Blade.Number/2*(1-r)./(r.*sin(phi))));

        % Aerodynamic force coefficients
        alpha = phi*180/pi - Control.Pitch.Fine - Blade.Twist;
        for j = 1:length(Blade.Radius)
            [~,ia] = unique(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(j))});
            Cl(j) = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(j))}(ia), Airfoil.Cl{Blade.IFoil(Blade.NFoil(j))}(ia), alpha(j));
            Cd(j) = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(j))}(ia), Airfoil.Cd{Blade.IFoil(Blade.NFoil(j))}(ia), alpha(j));
        end
        Cl = Cl(:);
        Cd = Cd(:);
        Cl(isnan(Cl)) = 1e-6;
        Cd(isnan(Cd)) = 0;
        Cl(Cl == 0) = 1e-6;

        % Thrust coefficient
        CT = sigmar.*(1-a).^2.*(Cl.*cos(phi)+Cd.*sin(phi))./(sin(phi).^2);

        % New induction factors
        for j = 1:length(Blade.Radius)
            if CT(j) < 0.96
                anew(j) = 1/(1+4*F(j)*sin(phi(j))^2/(sigmar(j)*Cl(j)*cos(phi(j))));
            else
                anew(j) = (1/F(j))*(0.143+sqrt(0.0203-0.6427*(0.889-CT(j))));
            end
        end
        a_new = 1./(4*F.*cos(phi)./(sigmar.*Cl)-1);

        if max(abs(anew - a)) < 0.01 && max(abs(a_new - a_)) < 0.002

            % Power coefficient
            dCP = (8*(lambdar-[0; lambdar(1:end-1)])./TSR(i)^2).*F.*sin(phi).^2.*(cos(phi)-lambdar.*sin(phi)).*(sin(phi)+lambdar.*cos(phi)).*(1-(Cd./Cl).*(cot(phi))).*lambdar.^2;
            CP(i) = real(sum(dCP(~isnan(dCP))));

            break
        end

    end

    if CP(i) - CP(i-1) < 0
        break
    end

end

% RPM range
TSRopt = TSRi(CP==max(CP));
TSRopt = TSRopt(1);
RPM = TSRopt*WindSpeeds/Blade.Radius(end) * 60/(2*pi);
RPM(RPM < Control.Torque.SpeedB/Drivetrain.Gearbox.Ratio) = Control.Torque.SpeedB/Drivetrain.Gearbox.Ratio;
RPM(RPM > Control.Torque.SpeedC/Drivetrain.Gearbox.Ratio) = Control.Torque.SpeedC/Drivetrain.Gearbox.Ratio;
RPM(WindSpeeds < Control.WindSpeed.Cutin) = 0;
RPM(WindSpeeds > Control.WindSpeed.Cutout) = 0;

% Power coefficients over the operating range
WindSpeeds = 0:0.1:(5*ceil(Control.WindSpeed.Cutout/5) + 5);
TSR = RPM*(2*pi/60)*Blade.Radius(end)./WindSpeeds;
CP = interp1(TSRi, CP, TSR);

% Rated power
Prated = Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Gearbox.Efficiency * Drivetrain.Generator.Efficiency;

% Find rated wind speed
Urated = min(WindSpeeds(P >= Prated));
CPmax = max(CP);
TSR_opt = TSR(CP == CPmax);

% Drive train
SpeedRange = [Control.Torque.SpeedB, Control.Torque.SpeedC]/Drivetrain.Gearbox.Ratio;
if Drivetrain.Gearbox.Ratio == 1
    DrivetrainType = 'Direct drive system';
else
    DrivetrainType = 'Geared system';
end

% Update tables
set(handles.Performance, 'Data', { ...
    'Rated power:', [num2str(Drivetrain.Generator.Power/1e6, '%2.1f'), ' MW']; ...
    'Cut-in wind speed:', [num2str(Control.WindSpeed.Cutin, '%2.1f'), ' m/s']; ...
    'Rated wind speed:', [num2str(Urated, '%2.1f'), ' m/s']; ...
    'Cut-out wind speed:', [num2str(Control.WindSpeed.Cutout, '%2.1f'), ' m/s']; ...
    'Speed range:', [num2str(SpeedRange(1), '%2.1f'), ' - ', num2str(SpeedRange(2), '%2.1f'), ' rpm']; ...
    'Optimal tip speed ratio:', num2str(TSR_opt, '%2.1f'); ...
    'Peak power coefficient:', num2str(CPmax, '%2.3f')});
set(handles.Configuration, 'Data', { ...
    'Number of blades:', int2str(Blade.Number); ...
    'Drivetrain:', DrivetrainType; ...
    'Control:', 'Variable speed, pitch-regulated'});
set(handles.Dimensions, 'Data', { ...
    'Rotor diameter:', [num2str(2*Blade.Radius(end), '%2.1f'), ' m']; ...
    'Hub height:', [num2str(Tower.HubHeight, '%2.1f'), ' m']; ...
    'Rotor mass:', [num2str(Blade.Number*trapz(Blade.Radius,Blade.Mass)/1e3+Nacelle.Hub.Mass/1e3, '%2.0f'), ' t (', num2str(trapz(Blade.Radius,Blade.Mass)/1e3, '%2.0f'), ' t per blade)']; ...
    'Nacelle mass:', [num2str(Nacelle.Housing.Mass/1e3, '%2.0f'), ' t']; ...
    'Tower mass:', [num2str(trapz(Tower.Height,Tower.Mass)/1e3, '%2.0f'), ' t']});

% Halt window
uiwait(handles.SummaryInfo);

%% Output function
function varargout = SummaryInfo_OutputFcn(hObject, eventdata, handles) 

% Close figure
delete(hObject)

%% OK button
function OK_Callback(hObject, eventdata, handles)
uiresume(handles.SummaryInfo);
