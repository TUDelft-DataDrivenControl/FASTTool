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

% Find rated wind speed by comparing demanded torque at rotational speed C
disp('Searching for rated wind speed...')
OmegaC = Control.Torque.SpeedC*2*pi/60;
U1 = Control.WindSpeed.Cutin;
U2 = Control.WindSpeed.Cutout;
[~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/U1);
Qr1 = 0.5*CQr*1.225*U1^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
[~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/U2);
Qr2 = 0.5*CQr*1.225*U2^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;

if Qr1 > Control.Torque.Demanded
    Urated = Control.WindSpeed.Cutin;
elseif Qr2 < Control.Torque.Demanded
    Urated = Control.WindSpeed.Cutout;
else
    success = false;
    for iter = 1:100
        U_new = U1 + abs((Control.Torque.Demanded-Qr1)/(Qr1-Qr2))*(U2-U1);
        [~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/U_new);
        Qr_new = 0.5*CQr*1.225*U_new^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
        if abs((Qr_new - Control.Torque.Demanded)/Control.Torque.Demanded) < 0.005
            success = true;
            break
        end
        if Qr_new < Control.Torque.Demanded
           U1 = U_new;
           Qr1 = Qr_new;
        else
           U2 = U_new;
           Qr2 = Qr_new;
        end
    end
    Urated = U_new;
    if ~success
        warning('Tolerance not met during iteration of axial induction factor')
    end
end

% Find maximum power coefficient and corresponding tip speed ratio
disp('Searching for maximum power coefficient and optimal tip speed ratio...')
TSR1 = 1;
TSR2 = 20;
[~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, TSR1);
CP1 = TSR1*CQr;
[~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, TSR2);
CP2 = TSR2*CQr;

success = false;
for iter = 1:100
    TSR_new = (TSR1 + TSR2)/2;
    [~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, TSR_new);
    CP_new = TSR_new*CQr;
    if abs((TSR2 - TSR_new)/TSR2) < 0.005
        success = true;
        break
    end
    if CP1 < CP2
       TSR1 = TSR_new;
       CP1 = CP_new;
    else
       TSR2 = TSR_new;
       CP2 = CP_new;
    end
end
CPmax = CP_new;
if ((TSR_new - 1) / 1.01) < 0.005
    OptimalTSR = '< 1';
elseif ((20 - TSR_new) / 20) < 0.005
    OptimalTSR = '> 20';
else
    OptimalTSR = num2str(TSR_new, '%2.1f');
end
if ~success
    warning('Tolerance not met during iteration of axial induction factor')
end

% Determine rated electrical power
Prated = Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Generator.Efficiency;

% Drive train
SpeedRange = [Control.Torque.SpeedB, Control.Torque.SpeedC]/Drivetrain.Gearbox.Ratio;
if Drivetrain.Gearbox.Ratio == 1
    DrivetrainType = 'Direct drive system';
else
    DrivetrainType = 'Geared system';
end

% Update tables
set(handles.Performance, 'Data', { ...
    'Rated power:', [num2str(Prated/1e6, '%2.2f'), ' MW']; ...
    'Cut-in wind speed:', [num2str(Control.WindSpeed.Cutin, '%2.1f'), ' m/s']; ...
    'Rated wind speed:', [num2str(Urated, '%2.1f'), ' m/s']; ...
    'Cut-out wind speed:', [num2str(Control.WindSpeed.Cutout, '%2.1f'), ' m/s']; ...
    'Speed range:', [num2str(SpeedRange(1), '%2.1f'), ' - ', num2str(SpeedRange(2), '%2.1f'), ' rpm']; ...
    'Optimal tip speed ratio:', OptimalTSR; ...
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
