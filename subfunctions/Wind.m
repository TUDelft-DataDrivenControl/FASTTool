%% Initialization code
function varargout = Wind(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Wind_OpeningFcn, ...
                   'gui_OutputFcn',  @Wind_OutputFcn, ...
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
function Wind_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.CertificationSettings = varargin{1};
handles.Blade = varargin{2};
handles.Tower = varargin{3};
handles.Input = varargin;

% Update input fields
set(handles.WindLength_textbox, 'String', num2str(handles.CertificationSettings.Wind.T));
set(handles.WindWidth_textbox, 'String', num2str(handles.CertificationSettings.Wind.Ly));
set(handles.WindHeight_textbox, 'String', num2str(handles.CertificationSettings.Wind.Lz));
set(handles.WindStep_textbox, 'String', num2str(handles.CertificationSettings.Wind.dt));
set(handles.WindNy_textbox, 'String', num2str(handles.CertificationSettings.Wind.Ny));
set(handles.WindNz_textbox, 'String', num2str(handles.CertificationSettings.Wind.Nz));
set(handles.WindMode, 'Value', handles.CertificationSettings.Wind.Type);
set(handles.WindClass, 'Value', handles.CertificationSettings.Wind.Class(1));
set(handles.WindTurbulence, 'Value', handles.CertificationSettings.Wind.Class(2));
if get(handles.WindMode, 'Value') == 1       % Steady wind
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(handles.WindMode, 'Value') == 2   % Stepped wind
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.Step))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Start at:')
elseif get(handles.WindMode, 'Value') == 3   % Normal wind profile (NWP)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(handles.WindMode, 'Value') == 4   % Normal turbulence model (NTM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(handles.WindMode, 'Value') == 5   % Annual extreme wind speed (EWM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(handles.WindMode, 'Value') == 6   % 50-year extreme wind speed (EWM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(handles.WindMode, 'Value') == 7   % Extreme wind shear (EWS)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.EWS))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
elseif get(handles.WindMode, 'Value') == 8   % Extreme turbulence model (ETM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(handles.WindMode, 'Value') == 9   % Extreme operating gust (EOG)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.EOG))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
elseif get(handles.WindMode, 'Value') == 10   % Extreme direction change (EDC)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.EDC))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
elseif get(handles.WindMode, 'Value') == 11   % Extreme coherent gust (ECG)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.ECG))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
end

% Update handles structure
guidata(hObject, handles);

% Draw
DrawGrid(handles)
DrawWind(handles)

% Halt window
uiwait(handles.Wind);

%% Closing function
function Wind_CloseRequestFcn(hObject, eventdata, handles)
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

%% Save button
function Apply_Callback(hObject, eventdata, handles)
handles.Save = true;
guidata(hObject, handles);
uiresume(handles.Wind);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.Wind);

%% Output function
function varargout = Wind_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save

    handles.CertificationSettings.Wind.Type = get(handles.WindMode, 'Value');
    handles.CertificationSettings.Wind.Class(1) = get(handles.WindClass, 'Value');
    handles.CertificationSettings.Wind.Class(2) = get(handles.WindTurbulence, 'Value');
    varargout{1} = handles.CertificationSettings;
    
else
    
    varargout = handles.Input;
    
end

% Close figure
delete(hObject)

%% Domain length - text box
function WindLength_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Wind.T))
end
handles.CertificationSettings.Wind.T = str2double(get(hObject,'String'));
guidata(hObject, handles);
DrawWind(handles)
function WindLength_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Domain width - text box
function WindWidth_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Wind.Ly))
end
handles.CertificationSettings.Wind.Ly = str2double(get(hObject,'String'));
guidata(hObject, handles);
DrawGrid(handles)
function WindWidth_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Domain height - text box
function WindHeight_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Wind.Lz))
end
handles.CertificationSettings.Wind.Lz = str2double(get(hObject,'String'));
guidata(hObject, handles);
DrawGrid(handles)
function WindHeight_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Time step - text box
function WindStep_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Wind.dt))
end
handles.CertificationSettings.Wind.dt = str2double(get(hObject,'String'));
guidata(hObject, handles);
DrawWind(handles)
function WindStep_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Number of elements in y-direction - text box
function WindNy_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Wind.Ny))
end
handles.CertificationSettings.Wind.Ny = ceil(str2double(get(hObject,'String')));
guidata(hObject, handles);
DrawGrid(handles)
function WindNy_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Number of elements in z-direction - text box
function WindNz_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Wind.Nz))
end
handles.CertificationSettings.Wind.Nz = ceil(str2double(get(hObject,'String')));
guidata(hObject, handles);
DrawGrid(handles)
function WindNz_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Draw grid
function DrawGrid(handles)

% Plot settings
axes(handles.GridPlot);
cla reset;
set(gca, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'Fontsize', 8);
axis equal
% xlim([-1 1]*50*ceil(max(get(gca,'XTick'))/50))
% ylim([0 1]*50*ceil(max(get(gca,'YTick'))/50))
xlabel('y [m]')
ylabel('z [m]')

% Plot grid
hold on
y = linspace(handles.CertificationSettings.Wind.Ly/2,-handles.CertificationSettings.Wind.Ly/2,handles.CertificationSettings.Wind.Ny);
z = handles.Tower.HubHeight + linspace(-handles.CertificationSettings.Wind.Lz/2,handles.CertificationSettings.Wind.Lz/2,handles.CertificationSettings.Wind.Nz);
for i = 2:handles.CertificationSettings.Wind.Ny-1
    plot([y(i) y(i)], [z(1) z(end)], 'Color', [0.75 0.75 0.75])
end
plot([y(1) y(1)], [z(1) z(end)], 'Color', [0.5 0.5 0.5])
plot([y(end) y(end)], [z(1) z(end)], 'Color', [0.5 0.5 0.5])
for i = 2:handles.CertificationSettings.Wind.Nz-1
    plot([y(1) y(end)], [z(i) z(i)], 'Color', [0.75 0.75 0.75])
end
plot([y(1) y(end)], [z(1) z(1)], 'Color', [0.5 0.5 0.5])
plot([y(1) y(end)], [z(end) z(end)], 'Color', [0.5 0.5 0.5])
        
% Plot tower
x = [-0.5*handles.Tower.Diameter; 0.5*flip(handles.Tower.Diameter)];
y = [handles.Tower.Height; flip(handles.Tower.Height)];
patch(x, y, [0 165, 213]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

% Plot rotor
plot(handles.Blade.Radius(end)*cos(linspace(0,2*pi,25)), handles.Tower.HubHeight + handles.Blade.Radius(end)*sin(linspace(0,2*pi,25)), '--', 'Color', [0 165, 213]/255)
x = [-handles.Blade.PitchAxis.*handles.Blade.Chord; flip((1-handles.Blade.PitchAxis).*handles.Blade.Chord)];
y = [handles.Blade.Radius; flip(handles.Blade.Radius)];
A = [x'; y'];
if handles.Blade.Number == 1 || handles.Blade.Number == 2 || handles.Blade.Number == 4
    R = [cos(pi/4),-sin(pi/4); ...
         sin(pi/4), cos(pi/4)];
    A = R*A;
end
t = -2*pi/handles.Blade.Number;
R = [cos(t),-sin(t); ...
     sin(t), cos(t)];
x = [];
y = [];
for i = 1:handles.Blade.Number   
    A = R*A;
    x = [x, A(1,:)];
    y = [y, A(2,:)];
end
patch(-x, handles.Tower.HubHeight + y, [0 165, 213]/255, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
hold off

%% Wind type
function WindMode_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1       % Steady wind
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 2   % Stepped wind
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.Step))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Start at:')
    set(handles.Event_unit, 'String', 'm/s')
elseif get(hObject, 'Value') == 3   % Normal wind profile (NWP)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 4   % Normal turbulence model (NTM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 5   % Annual extreme wind speed (EWM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 6   % 50-year extreme wind speed (EWM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 7   % Extreme wind shear (EWS)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.EWS))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
    set(handles.Event_unit, 'String', 's')
elseif get(hObject, 'Value') == 8   % Extreme turbulence model (ETM)
    set(handles.WindEvent, 'Visible', 'off')
    set(handles.Event_label, 'Visible', 'off')
    set(handles.Event_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 9   % Extreme operating gust (EOG)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.EOG))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
    set(handles.Event_unit, 'String', 's')
elseif get(hObject, 'Value') == 10   % Extreme direction change (EDC)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.EDC))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
    set(handles.Event_unit, 'String', 's')
elseif get(hObject, 'Value') == 11   % Extreme coherent gust (ECG)
    set(handles.WindEvent, 'Visible', 'on')
    set(handles.WindEvent, 'String', num2str(handles.CertificationSettings.Wind.ECG))
    set(handles.Event_label, 'Visible', 'on')
    set(handles.Event_unit, 'Visible', 'on')
    set(handles.Event_label, 'String', 'Event at:')
    set(handles.Event_unit, 'String', 's')
end
DrawWind(handles)
function WindMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% IEC class 1/2/3
function WindClass_Callback(hObject, eventdata, handles)
DrawWind(handles)
function WindClass_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% IEC class A/B/C
function WindTurbulence_Callback(hObject, eventdata, handles)
DrawWind(handles)
function WindTurbulence_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Event time stamp - text box
function WindEvent_Callback(hObject, eventdata, handles)
if get(handles.WindMode, 'Value') == 2   % Stepped wind
    event = handles.CertificationSettings.Wind.Step;
elseif get(handles.WindMode, 'Value') == 7   % Extreme wind shear (EWS)
    event = handles.CertificationSettings.Wind.EWS;
elseif get(handles.WindMode, 'Value') == 9   % Extreme operating gust (EOG)
    event = handles.CertificationSettings.Wind.EOG;
elseif get(handles.WindMode, 'Value') == 10   % Extreme direction change (EDC)
    event = handles.CertificationSettings.Wind.EDC;
elseif get(handles.WindMode, 'Value') == 11   % Extreme coherent gust (ECG)
    event = handles.CertificationSettings.Wind.ECG;
end
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(event))
end
if get(handles.WindMode, 'Value') == 2   % Stepped wind
    handles.CertificationSettings.Wind.Step = str2double(get(hObject,'String'));
elseif get(handles.WindMode, 'Value') == 7   % Extreme wind shear (EWS)
    handles.CertificationSettings.Wind.EWS = str2double(get(hObject,'String'));
elseif get(handles.WindMode, 'Value') == 9   % Extreme operating gust (EOG)
    handles.CertificationSettings.Wind.EOG = str2double(get(hObject,'String'));
elseif get(handles.WindMode, 'Value') == 10   % Extreme direction change (EDC)
    handles.CertificationSettings.Wind.EDC = str2double(get(hObject,'String'));
elseif get(handles.WindMode, 'Value') == 11   % Extreme coherent gust (ECG)
    handles.CertificationSettings.Wind.ECG = str2double(get(hObject,'String'));
end
guidata(hObject, handles);
DrawWind(handles)
function WindEvent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Draw wind
function DrawWind(handles)

% Plot settings
axes(handles.WindPlot);
cla reset;

% Time series
U = handles.CertificationSettings.Run.WindSpeed(1);
t = 0:handles.CertificationSettings.Wind.dt:handles.CertificationSettings.Wind.T;

% IEC class
if get(handles.WindTurbulence, 'Value') == 1
    Iref = 0.16;
elseif get(handles.WindTurbulence, 'Value') == 2
    Iref = 0.14;
elseif get(handles.WindTurbulence, 'Value') == 3
    Iref = 0.12;
end
if get(handles.WindClass, 'Value') == 1
    Vref = 50;
elseif get(handles.WindClass, 'Value') == 2
    Vref = 42.5;
elseif get(handles.WindClass, 'Value') == 3
    Vref = 37.5;
end
if handles.Tower.HubHeight >= 60
    Lambda1 = 0.7*handles.Tower.HubHeight;
else
    Lambda1 = 42;
end

% Plot wind
if get(handles.WindMode, 'Value') == 1       % Steady wind
    
    % Constant wind
    u = 0*t + U;
    plot(t,u);

elseif get(handles.WindMode, 'Value') == 2   % Stepped wind
    
    % Stepped profile
    u = handles.CertificationSettings.Wind.Step-1+ceil((U-handles.CertificationSettings.Wind.Step+1)*t/handles.CertificationSettings.Wind.T);
    u(1) = handles.CertificationSettings.Wind.Step;
    plot(t,u);
    
elseif get(handles.WindMode, 'Value') == 3   % Normal wind profile (NWP)
    
    % Constant wind
    u1 = 0*t + U*((handles.Tower.HubHeight+handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2;
    u2 = 0*t + U;
    u3 = 0*t + U*((handles.Tower.HubHeight-handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2;
    
    % Plot
    plot(t,u1,t,u2,t,u3);
    legend({'u(z_{hub} + R)','u(z_{hub})','u(z_{hub} - R)'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 4   % Normal turbulence model (NTM)
    
    % NTM values
    sigma1 = Iref*(0.75*U+5.6);
    sigma2 = 0.8*sigma1;
    sigma3 = 0.5*sigma1;
    Lu = 8.1*Lambda1;
    Lv = 2.7*Lambda1;
    Lw = 0.66*Lambda1;
    
    % Wave number vector
    L = U * handles.CertificationSettings.Wind.T;
    N = length(t);
    m = ifftshift(-N/2:N/2-1);
    k = 2*pi*m/L;
    
    % Spectrum
    Fu = sigma1^2 * 4*Lu/U./(1+6/(2*pi)*abs(k)*Lu).^(5/3);
    Fv = sigma2^2 * 4*Lv/U./(1+6/(2*pi)*abs(k)*Lv).^(5/3);
    Fw = sigma3^2 * 4*Lw/U./(1+6/(2*pi)*abs(k)*Lw).^(5/3);
    n = randn([3,length(k)]) + sqrt(-1)*randn([3,length(k)]);
    dZ = sqrt(2*pi*[Fu;Fv;Fw]/L) .* n;
    
    % IFFT
    u = N*real(ifft(dZ(1,:)));
    v = N*real(ifft(dZ(2,:)));
    w = N*real(ifft(dZ(3,:)));
    u = (u -mean(u)) * sigma1/std(u) + U;
    v = (v -mean(v)) * sigma2/std(v);
    w = (w -mean(w)) * sigma3/std(w);
    
    % Plot
    plot(t,u,t,v,t,w);
    legend({'u','v','w'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 5   % Annual extreme wind speed (EWM)
    
    % EWM values
    V50 = Vref;
    V1 = 0.8*V50;
    sigma1 = 0.11*V1;
    sigma2 = 0.8*sigma1;
    sigma3 = 0.5*sigma1;
    Lu = 8.1*Lambda1;
    Lv = 2.7*Lambda1;
    Lw = 0.66*Lambda1;
    
    % Wave number vector
    L = U * handles.CertificationSettings.Wind.T;
    N = length(t);
    m = ifftshift(-N/2:N/2-1);
    k = 2*pi*m/L;
    
    % Spectrum
    Fu = sigma1^2 * 4*Lu/U./(1+6/(2*pi)*abs(k)*Lu).^(5/3);
    Fv = sigma2^2 * 4*Lv/U./(1+6/(2*pi)*abs(k)*Lv).^(5/3);
    Fw = sigma3^2 * 4*Lw/U./(1+6/(2*pi)*abs(k)*Lw).^(5/3);
    n = randn([3,length(k)]) + sqrt(-1)*randn([3,length(k)]);
    dZ = sqrt(2*pi*[Fu;Fv;Fw]/L) .* n;
    
    % IFFT
    u = N*real(ifft(dZ(1,:)));
    v = N*real(ifft(dZ(2,:)));
    w = N*real(ifft(dZ(3,:)));
    u = (u -mean(u)) * sigma1/std(u) + V1;
    v = (v -mean(v)) * sigma2/std(v);
    w = (w -mean(w)) * sigma3/std(w);
    
    % Plot
    plot(t,u,t,v,t,w);
    legend({'u','v','w'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 6   % 50-year extreme wind speed (EWM)
    
    % EWM values
    V50 = Vref;
    sigma1 = 0.11*V50;
    sigma2 = 0.8*sigma1;
    sigma3 = 0.5*sigma1;
    Lu = 8.1*Lambda1;
    Lv = 2.7*Lambda1;
    Lw = 0.66*Lambda1;
    
    % Wave number vector
    L = U * handles.CertificationSettings.Wind.T;
    N = length(t);
    m = ifftshift(-N/2:N/2-1);
    k = 2*pi*m/L;
    
    % Spectrum
    Fu = sigma1^2 * 4*Lu/U./(1+6/(2*pi)*abs(k)*Lu).^(5/3);
    Fv = sigma2^2 * 4*Lv/U./(1+6/(2*pi)*abs(k)*Lv).^(5/3);
    Fw = sigma3^2 * 4*Lw/U./(1+6/(2*pi)*abs(k)*Lw).^(5/3);
    n = randn([3,length(k)]) + sqrt(-1)*randn([3,length(k)]);
    dZ = sqrt(2*pi*[Fu;Fv;Fw]/L) .* n;
    
    % IFFT
    u = N*real(ifft(dZ(1,:)));
    v = N*real(ifft(dZ(2,:)));
    w = N*real(ifft(dZ(3,:)));
    u = (u -mean(u)) * sigma1/std(u) + V50;
    v = (v -mean(v)) * sigma2/std(v);
    w = (w -mean(w)) * sigma3/std(w);
    
    % Plot
    plot(t,u,t,v,t,w);
    legend({'u','v','w'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 7   % Extreme wind shear (EWS)
       
    % NTM values
    sigma1 = Iref*(0.75*U+5.6);
    
    % Wind speeds
    beta = 6.4;
    Tg = 12;
    t0 = handles.CertificationSettings.Wind.EWS;
    u1 = U * ((handles.Tower.HubHeight+handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2 + 0.5*(2.5+0.2*beta*sigma1*(2*handles.Blade.Radius(end)/Lambda1)^0.25) * (1-cos(2*pi*(t-t0)/Tg));
    u2 = 0*t + U;
    u3 = U * ((handles.Tower.HubHeight-handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2 - 0.5*(2.5+0.2*beta*sigma1*(2*handles.Blade.Radius(end)/Lambda1)^0.25) * (1-cos(2*pi*(t-t0)/Tg));
    u1(t < t0) = U * ((handles.Tower.HubHeight+handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2;
    u1(t > t0+Tg) = U * ((handles.Tower.HubHeight+handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2;
    u3(t < t0) = U * ((handles.Tower.HubHeight-handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2;
    u3(t > t0+Tg) = U * ((handles.Tower.HubHeight-handles.Blade.Radius(end))/handles.Tower.HubHeight).^0.2;
    
    % Plot
    plot(t,u1,t,u2,t,u3);
    legend({'u(z_{hub} + R)','u(z_{hub})','u(z_{hub} - R)'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 8   % Extreme turbulence model (ETM)
    
    % ETM values
    Vave = 0.2*Vref;
    c = 2;
    sigma1 = c * Iref*(0.072*(Vave/c+3)*(U/c-4) + 10);
    sigma2 = 0.8*sigma1;
    sigma3 = 0.5*sigma1;
    Lu = 8.1*Lambda1;
    Lv = 2.7*Lambda1;
    Lw = 0.66*Lambda1;
    
    % Wave number vector
    L = U * handles.CertificationSettings.Wind.T;
    N = length(t);
    m = ifftshift(-N/2:N/2-1);
    k = 2*pi*m/L;
    
    % Spectrum
    Fu = sigma1^2 * 4*Lu/U./(1+6/(2*pi)*abs(k)*Lu).^(5/3);
    Fv = sigma2^2 * 4*Lv/U./(1+6/(2*pi)*abs(k)*Lv).^(5/3);
    Fw = sigma3^2 * 4*Lw/U./(1+6/(2*pi)*abs(k)*Lw).^(5/3);
    n = randn([3,length(k)]) + sqrt(-1)*randn([3,length(k)]);
    dZ = sqrt(2*pi*[Fu;Fv;Fw]/L) .* n;
    
    % IFFT
    u = N*real(ifft(dZ(1,:)));
    v = N*real(ifft(dZ(2,:)));
    w = N*real(ifft(dZ(3,:)));
    u = (u -mean(u)) * sigma1/std(u) + U;
    v = (v -mean(v)) * sigma2/std(v);
    w = (w -mean(w)) * sigma3/std(w);
    
    % Plot
    plot(t,u,t,v,t,w);
    legend({'u','v','w'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 9   % Extreme operating gust (EOG)
    
    % EWM values
    Ve50 = 1.4*Vref;
    Ve1 = 0.8*Ve50;
    
    % NTM values
    sigma1 = Iref*(0.75*U+5.6);
    
    % Gust at t0
    Vgust = min([1.35*(Ve1-U), 3.3*sigma1/(1+0.1*2*handles.Blade.Radius(end)/Lambda1)]);
    t0 = handles.CertificationSettings.Wind.EOG;
    Tg = 10.5;
    u = 0*t + U - 0.37*Vgust*sin(3*pi*(t-t0)/Tg) .* (1-cos(2*pi*(t-t0)/Tg));
    u(t < t0) = U;
    u(t > t0+Tg) = U;
    
    % Plot
    plot(t,u);
    
elseif get(handles.WindMode, 'Value') == 10   % Extreme direction change (EDC)
    
    % NTM values
    sigma1 = Iref*(0.75*U+5.6);
    
    % EDC values
    t0 = handles.CertificationSettings.Wind.EDC;
    Tg = 6;
    thetae = 4*atan(sigma1/(U*(1+0.1*2*handles.Blade.Radius(end)/Lambda1)));
    theta = 0.5*thetae * (1-cos(pi*(t-t0)/Tg));
    theta(t < t0) = 0;
    theta(t > t0+Tg) = thetae;
    
    % Wind speeds
    u = U*cos(theta);
    v = U*sin(theta);
    w = 0*t;
    
    % Plot
    plot(t,u,t,v,t,w);
    legend({'u','v','w'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
elseif get(handles.WindMode, 'Value') == 11   % Extreme coherent gust (ECG)
    
    % ECG values
    t0 = handles.CertificationSettings.Wind.ECG;
    Tg = 10;
    Vcg = 15;
    V = U + 0.5*Vcg*(1-cos(pi*(t-t0)/Tg));
    V(t < t0) = U;
    V(t > t0+Tg) = U + Vcg;
    if U < 4
        thetacg = pi;
    else
        thetacg = 4*pi/U;
    end
    theta = 0.5*thetacg*(1-cos(pi*(t-t0)/Tg));
    theta(t < t0) = 0;
    theta(t > t0+Tg) = thetacg;
        
    % Wind speeds
    u = V.*cos(theta);
    v = V.*sin(theta);
    w = 0*t;
    
    % Plot
    plot(t,u,t,v,t,w);
    legend({'u','v','w'}, 'Location', 'NorthEast', 'Color', 'none', 'Box', 'off')
    
end

% Update axes
set(gca, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'Fontsize', 8);
xlabel('t [s]')
ylabel('u [m/s]')
