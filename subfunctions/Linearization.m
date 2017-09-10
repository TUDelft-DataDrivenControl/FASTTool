%% Initialization code
function varargout = Linearization(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Linearization_OpeningFcn, ...
                   'gui_OutputFcn',  @Linearization_OutputFcn, ...
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
function Linearization_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Blade = varargin{1};
handles.Airfoil = varargin{2};
handles.Tower = varargin{3};
handles.Nacelle = varargin{4};
handles.Control = varargin{5};
handles.Drivetrain = varargin{6};
handles.LinModel = varargin{7};
handles.Input = varargin{7};

% Check if current linearized model is valid
if handles.LinModel
    if exist(handles.LinModel, 'file') == 2

        contents = whos('-file', handles.LinModel);
        if ismember('Lin', {contents.name}) && ismember('sysm', {contents.name})
            set(handles.FileCheck, 'String', 'Valid model found.')
        else
            set(handles.FileCheck, 'String', 'Invalid model.')
        end

    else

        set(handles.FileCheck, 'String', '<empty>')

    end
else
       
    set(handles.FileCheck, 'String', '<empty>')
    
end

% Propose values for the wind speed range
set(handles.WindSpeed_From, 'String', '5')
set(handles.WindSpeed_To, 'String', '25')
set(handles.WindSpeed_Step, 'String', '1')
set(handles.LinAmount, 'String', '1')
set(handles.LinRotations, 'String', '1')

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.Linearization);

%% Closing function
function Linearization_CloseRequestFcn(hObject, eventdata, handles)
button = questdlg('Save changes?');
if strcmp(button, 'Yes')
    handles.Save = true;
    guidata(hObject, handles);
    uiresume(hObject);
elseif strcmp(button, 'No')
    handles.Save = false;
    guidata(hObject, handles);
    uiresume(hObject);
end

%% Apply button
function Apply_Callback(hObject, eventdata, handles)
handles.Save = true;
guidata(hObject, handles);
uiresume(handles.Linearization);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.Linearization);

%% Output function
function varargout = Linearization_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save
    varargout{1} = handles.LinModel;
else
    varargout{1} = handles.Input;
end

% Apply figure
delete(hObject)

%% Set file names
function SetLinModel_Callback(hObject, eventdata, handles)

% Set file name
[FileName,PathName] = uiputfile('*.mat', 'Set file name for linearized model');
handles.LinModel = [PathName, FileName];

if FileName

    % Enable buttons
    set(handles.Linearize, 'Enable', 'on');
    
    % Store in handles
    guidata(hObject, handles);

end

%% Open existing model
function OpenLinModel_Callback(hObject, eventdata, handles)

% Get file name
[FileName,PathName] = uigetfile('*.mat', 'Select existing linearized model');
handles.LinModel = [PathName, FileName];

% Check if it is a valid model
if FileName
    if exist(handles.LinModel, 'file') == 2

        contents = whos('-file', handles.LinModel);
        if ismember('Lin', {contents.name}) && ismember('sysm', {contents.name})
            set(handles.FileCheck, 'String', 'Valid model found.')
        else
            set(handles.FileCheck, 'String', 'Invalid model.')
        end
    end
else

    set(handles.FileCheck, 'String', '<empty>')

end

% Update handles
guidata(hObject, handles);
    
%% Linearization
function Linearize_Callback(hObject, eventdata, handles)

% Wind speed range
WindSpeeds = str2double(get(handles.WindSpeed_From, 'String')):str2double(get(handles.WindSpeed_Step, 'String')):str2double(get(handles.WindSpeed_To, 'String'));
if sum(WindSpeeds) == 0 || isnan(sum(WindSpeeds))
    errordlg('Invalid wind speed range.', 'Error')
else
LinAmount = str2double(get(handles.LinAmount, 'String'));
LinRotations = str2double(get(handles.LinRotations, 'String'));
set(hObject, 'Enable', 'off');

disp('Starting linearization...')
disp(' ')

% Get geometry from handles
Blade = handles.Blade;
Airfoil = handles.Airfoil;
Tower = handles.Tower;
Nacelle = handles.Nacelle;
Control = handles.Control;
Drivetrain = handles.Drivetrain;

% Run modal analysis for 0 rpm
disp('Finding non-rotating mode shapes...')

% Run BModes for tower
evalc([...
'[y11_shape, y11_coeff, y11_freq,', ...
' y12_shape, y12_coeff, y12_freq,', ...
' y21_shape, y21_coeff, y21_freq,', ...
' y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,2,0);']);

% Store in handles
Tower.ForeAft1_coeff = y21_coeff;
Tower.ForeAft2_coeff = y22_coeff;
Tower.SideSide1_coeff = y11_coeff;
Tower.SideSide2_coeff = y12_coeff;

% Run BModes for blade
evalc([...
'[y11_shape, y11_coeff, y11_freq,', ...
' y12_shape, y12_coeff, y12_freq,', ...
' y21_shape, y21_coeff, y21_freq,', ...
' y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,1,0);']);

% Store in handles
Blade.Flap1_coeff = y11_coeff;
Blade.Flap2_coeff = y12_coeff;
Blade.Edge1_coeff = y21_coeff;
Blade.Edge2_coeff = y22_coeff;

% Run steady BEM to find RPM range
disp('Looking for operating rotor speeds and trim cases...')

% Cp-lambda curve
TSRi = 0.1:0.1:20;
CP = zeros(size(TSRi));
for i = 2:length(TSRi)

    % Local tip speed ratio and solidity
    r = Blade.Radius/Blade.Radius(end);
    lambdar = TSRi(i) * r;
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
        alpha = phi*180/pi - Blade.Twist - Control.Pitch.Fine;
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

        % Thrust per rotor annulus
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
            dCP = (8*(lambdar-[0; lambdar(1:end-1)])./TSRi(i)^2).*F.*sin(phi).^2.*(cos(phi)-lambdar.*sin(phi)).*(sin(phi)+lambdar.*cos(phi)).*(1-(Cd./Cl).*(cot(phi))).*lambdar.^2;
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
TSR = RPM*(2*pi/60)*Blade.Radius(end)./WindSpeeds;
CP_U = interp1(TSRi, CP, TSR);

% Rated power
%Prated = Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Gearbox.Efficiency * Drivetrain.Generator.Efficiency;
Prated = Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Generator.Efficiency;

% Preliminary power curve
P = 0.5 * 1.225 * pi*Blade.Radius(end)^2 * WindSpeeds.^3 .* CP_U .* Drivetrain.Gearbox.Efficiency .* Drivetrain.Generator.Efficiency;

% Pitch range
disp('Finding pitch angles...')
PitchAngle = Control.Pitch.Fine*ones(size(WindSpeeds));
for u = 1:length(WindSpeeds)
    
    if P(u) >= Prated
    
    if u > 1
        Pitch = 0.1*floor(10*PitchAngle(u-1)):0.1:45;
    else
        Pitch = 0:0.1:45;
    end
    Pi = nan(size(Pitch));
    for i = 1:length(Pitch)

        % Local tip speed ratio and solidity
        r = Blade.Radius/Blade.Radius(end);
        lambdar = TSR(u) * r;
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
            alpha = phi*180/pi - Pitch(i) - Blade.Twist - Control.Pitch.Fine;
            for k = 1:length(Blade.Radius)
                [~,ia] = unique(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))});
                Cl(k) = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}(ia), Airfoil.Cl{Blade.IFoil(Blade.NFoil(k))}(ia), alpha(k));
                Cd(k) = interp1(Airfoil.Alpha{Blade.IFoil(Blade.NFoil(k))}(ia), Airfoil.Cd{Blade.IFoil(Blade.NFoil(k))}(ia), alpha(k));
            end
            Cl = Cl(:);
            Cd = Cd(:);
            Cl(isnan(Cl)) = 1e-6;
            Cd(isnan(Cd)) = 0;
            Cl(Cl == 0) = 1e-6;

            % Thrust per rotor annulus
            CT = sigmar.*(1-a).^2.*(Cl.*cos(phi)+Cd.*sin(phi))./(sin(phi).^2);

            % New induction factors
            for k = 1:length(Blade.Radius)
                if CT(k) < 0.96
                    anew(k) = 1/(1+4*F(k)*sin(phi(k))^2/(sigmar(k)*Cl(k)*cos(phi(k))));
                else
                    anew(k) = (1/F(k))*(0.143+sqrt(0.0203-0.6427*(0.889-CT(k))));
                end
            end
            a_new = 1./(4*F.*cos(phi)./(sigmar.*Cl)-1);

            if max(abs(anew - a)) < 0.01 && max(abs(a_new - a_)) < 0.002

                % Power coefficient
                dCP = (8*(lambdar-[0; lambdar(1:end-1)])./TSR(u)^2).*F.*sin(phi).^2.*(cos(phi)-lambdar.*sin(phi)).*(sin(phi)+lambdar.*cos(phi)).*(1-(Cd./Cl).*(cot(phi))).*lambdar.^2;
                CP(i) = real(sum(dCP(~isnan(dCP))));

                break
            elseif iter == 100
                warning('Tolerance not met during iterations')
            end

        end

        Pi(i) = 0.5 * 1.225 * pi*Blade.Radius(end)^2 * CP(i) * Drivetrain.Gearbox.Efficiency * Drivetrain.Generator.Efficiency * WindSpeeds(u)^3;
        
        if Pi(i) < Prated
            break
        end

    end
    
    Pitch = Pitch(~isnan(Pi));
    CP = CP(~isnan(Pi));
    Pi = Pi(~isnan(Pi));
    if length(Pitch) > 1
        PitchAngle(u) = interp1(Pi,Pitch,Prated);
        CP_U(u) = interp1(Pi,CP,Prated);
    else
        PitchAngle(u) = PitchAngle(u-1) + (PitchAngle(u-1)-PitchAngle(u-2));
        CP_U(u) = CP_U(u-1) + (CP_U(u-1)-CP_U(u-2));
    end
    
    end
    
end

% Avoid some common errors with linearization
% (1) Set the gearbox efficiency (to avoid error that ADAMS cannot handle
% nonideal gearboxes)
Drivetrain.Gearbox.Efficiency = 1;

% (2) Scale the HSS inertia from the reference machine (to avoid issues
% with reaching very high rotor speeds during the linearization of some
% direct-drive machines). Assuming that P ~ Mass, Radius² ~ P^(2/3), we get
% a relation in the shape of I/Iref = (P/Pref)^(5/3) * (RPM_ref/RPM)^2.
Drivetrain.Generator.HSSInertia = 534.116 * (Prated/(5e6))^(5/3) * (12.1*97/Control.Torque.SpeedC)^2;

% Turbine input file
disp('Writing input files...')

% Simulink settings
TSim = 60;
FAST_InputFileName = [pwd, '\subfunctions\inputfiles\FAST.fst'];

% Turbine input files
AeroDyn(Blade,Airfoil,Tower,'Linearize');
ServoDyn(Drivetrain,Control,'Linearize');

% Preload the OutList
load([pwd '\subfunctions\OutList.mat'])
assignin('base', 'OutList', OutList);

% Run linearization
sysm = cell(length(WindSpeeds),1);
for j = 1:length(WindSpeeds)
        
    % Status update
    disp(['Linearizing at U = ', num2str(WindSpeeds(j), '%5.2f'), ' m/s, ', num2str(RPM(j), '%5.2f'), ' rpm, ', num2str(PitchAngle(j), '%5.2f'), ' deg pitch'])
    
    % Set initial RPM and pitch angle in ElastoDyn input file
    ElastoDyn(Blade,Tower,Nacelle,Drivetrain,Control,'Linearize',RPM(j),PitchAngle(j));
    
    % Set linearization times for 10 deg azimuth step (after 30 s)
    LinAziPositions = linspace(0,360*LinRotations,LinAmount+1);
    LinTimes = TSim + Control.DT * round(LinAziPositions(2:end)/(RPM(j)*6) / Control.DT);
    TMax = max(LinTimes);
    
    FASTinput(Control.DT, TMax, 'Linearize', LinTimes);

    % Wind input file
    Wind.Type = 1;
    InflowWind(Wind,WindSpeeds(j),Tower.HubHeight,Blade.Radius(end))

    % Run FAST and prevent console output
    assignin('base', 'TMax', TMax);
    assignin('base', 'FAST_InputFileName', FAST_InputFileName);
    evalc('sim(''OpenLoop'',TMax);');

% Extract steady state solution after 60 seconds and average over 36
% steps
    A = 0;
    B = 0;
    C = 0;
    D = 0;
    for i = 1:LinAmount
        LinName = [pwd, '\subfunctions\inputfiles\FAST.SFunc.', int2str(i), '.lin'];
        data = ReadFASTLinear(LinName);
        A = A + 1/LinAmount * data.A;
        B = B + 1/LinAmount * data.B;
        C = C + 1/LinAmount * data.C;
        D = D + 1/LinAmount * data.D;
    end
    sysm{j} = ss(A,B,C,D,'InputName',data.u_desc,'Outputname',data.y_desc);
    
    Lin.V(j) = data.y_op{1};
    Lin.Torque(j) = data.y_op{5};
    Lin.Pitch(j) =  data.y_op{34}*pi/180;
    Lin.GSpeed(j) = data.y_op{33}*pi/30;
    Lin.RSpeed(j) = data.y_op{38}*pi/30;
end

save(handles.LinModel, 'sysm', 'Lin');
disp(' ')
disp(' ')
disp('... Linearization complete!')

set(hObject, 'Enable', 'on');

end

%% Wind speed steps
function WindSpeed_From_Callback(hObject, eventdata, handles)
function WindSpeed_From_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function WindSpeed_To_Callback(hObject, eventdata, handles)
function WindSpeed_To_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function WindSpeed_Step_Callback(hObject, eventdata, handles)
function WindSpeed_Step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LinAmount_Callback(hObject, eventdata, handles)
function LinAmount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LinRotations_Callback(hObject, eventdata, handles)
function LinRotations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
