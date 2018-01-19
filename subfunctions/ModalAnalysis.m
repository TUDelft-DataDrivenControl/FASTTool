%% Initialization code
function varargout = ModalAnalysis(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ModalAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @ModalAnalysis_OutputFcn, ...
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
function ModalAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Blade = varargin{1};
handles.Airfoil = varargin{2};
handles.Tower = varargin{3};
handles.Nacelle = varargin{4};
handles.Drivetrain = varargin{5};
handles.Control = varargin{6};

% Blank Campbell diagram
RPMend = handles.Control.Torque.SpeedC/handles.Drivetrain.Gearbox.Ratio;
RPMend = ceil(RPMend/5)*5;
RPMstep = str2double(get(handles.RPMstep_textbox, 'String'));
handles.RPMs = 0:RPMstep:RPMend;
handles.Tower_ForeAft1_freq = nan(size(handles.RPMs));
handles.Tower_ForeAft2_freq = nan(size(handles.RPMs));
handles.Tower_SideSide1_freq = nan(size(handles.RPMs));
handles.Tower_SideSide2_freq = nan(size(handles.RPMs));
handles.Blade_Flap1_freq = nan(size(handles.RPMs));
handles.Blade_Flap2_freq = nan(size(handles.RPMs));
handles.Blade_Edge1_freq = nan(size(handles.RPMs));
handles.Blade_Edge2_freq = nan(size(handles.RPMs));
DrawCampbell(handles,false);

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.ModalAnalysis);

%% Output function
function varargout = ModalAnalysis_OutputFcn(hObject, eventdata, handles) 

% Close figure
delete(hObject)

%% Close button
function Close_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume(handles.ModalAnalysis);

%% Modal analysis (BModes)
function RunModal_Callback(hObject, eventdata, handles)

% Disable window
buttons = findall(handles.ModalAnalysis, 'Type', 'UIControl');
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'off');
end
    
% Get geometry from handles
Blade = handles.Blade;
Tower = handles.Tower;
Nacelle = handles.Nacelle;
Control = handles.Control;

% Empty blade mode frequency vectors
handles.Blade_Flap1_freq = nan(size(handles.RPMs));
handles.Blade_Flap2_freq = nan(size(handles.RPMs));
handles.Blade_Edge1_freq = nan(size(handles.RPMs));
handles.Blade_Edge2_freq = nan(size(handles.RPMs));

% Run BModes for tower
[y11_shape, y11_coeff, y11_freq, ...
    y12_shape, y12_coeff, y12_freq, ...
    y21_shape, y21_coeff, y21_freq, ...
    y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,2,0);

% Store in handles
handles.Tower_ForeAft1_shape = y21_shape;
handles.Tower_ForeAft2_shape = y22_shape;
handles.Tower_ForeAft1_freq = y21_freq * ones(size(handles.RPMs));
handles.Tower_ForeAft2_freq = y22_freq * ones(size(handles.RPMs));
handles.Tower_ForeAft1_coeff = y21_coeff;
handles.Tower_ForeAft2_coeff = y22_coeff;
handles.Tower_SideSide1_shape = y11_shape;
handles.Tower_SideSide2_shape = y12_shape;
handles.Tower_SideSide1_freq = y11_freq * ones(size(handles.RPMs));
handles.Tower_SideSide2_freq = y12_freq * ones(size(handles.RPMs));
handles.Tower_SideSide1_coeff = y11_coeff;
handles.Tower_SideSide2_coeff = y12_coeff;

% Update Campbell diagram
DrawCampbell(handles,false)

% Cycle through RPM range
for rpm = 1:length(handles.RPMs)
    
    % Run BModes for blade
    [y11_shape, y11_coeff, y11_freq, ...
        y12_shape, y12_coeff, y12_freq, ...
        y21_shape, y21_coeff, y21_freq, ...
        y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,1,handles.RPMs(rpm));
    
    % Store mode frequencies in handles
    handles.Blade_Flap1_freq(rpm) = y11_freq;
    handles.Blade_Flap2_freq(rpm) = y12_freq;
    handles.Blade_Edge1_freq(rpm) = y21_freq;
    handles.Blade_Edge2_freq(rpm) = y22_freq;

    % Store 0 rpm results separately
    if handles.RPMs(rpm) == 0

        % Update fields
        set(handles.Tower_ForeAft1_textbox, 'String', num2str(handles.Tower_ForeAft1_freq(1), '%0.2f'));
        set(handles.Tower_ForeAft2_textbox, 'String', num2str(handles.Tower_ForeAft2_freq(1), '%0.2f'));
        set(handles.Tower_SideSide1_textbox, 'String', num2str(handles.Tower_SideSide1_freq(1), '%0.2f'));
        set(handles.Tower_SideSide2_textbox, 'String', num2str(handles.Tower_SideSide2_freq(1), '%0.2f'));
        set(handles.Blade_Flap1_textbox, 'String', num2str(handles.Blade_Flap1_freq(1), '%0.2f'));
        set(handles.Blade_Flap2_textbox, 'String', num2str(handles.Blade_Flap2_freq(1), '%0.2f'));
        set(handles.Blade_Edge1_textbox, 'String', num2str(handles.Blade_Edge1_freq(1), '%0.2f'));
        set(handles.Blade_Edge2_textbox, 'String', num2str(handles.Blade_Edge2_freq(1), '%0.2f'));

        % Store in handles
        handles.Blade_Flap1_shape = y11_shape;
        handles.Blade_Flap2_shape = y12_shape;
        handles.Blade_Edge1_shape = y21_shape;
        handles.Blade_Edge2_shape = y22_shape;
        handles.Blade_Flap1_coeff = y11_coeff;
        handles.Blade_Flap2_coeff = y12_coeff;
        handles.Blade_Edge1_coeff = y21_coeff;
        handles.Blade_Edge2_coeff = y22_coeff;

    end

    % Update Campbell diagram
    DrawCampbell(handles,false)

end

assignin('base', 'RotorSpeed', handles.RPMs(:));
assignin('base', 'TowerForeAft1', handles.Tower_ForeAft1_freq(:));
assignin('base', 'TowerForeAft2', handles.Tower_ForeAft2_freq(:));
assignin('base', 'TowerSideSide1', handles.Tower_SideSide1_freq(:));
assignin('base', 'TowerSideSide2', handles.Tower_SideSide2_freq(:));
assignin('base', 'BladeFlapwise1', handles.Blade_Flap1_freq(:));
assignin('base', 'BladeFlapwise2', handles.Blade_Flap2_freq(:));
assignin('base', 'BladeEdgewise1', handles.Blade_Edge1_freq(:));
assignin('base', 'BladeEdgewise2', handles.Blade_Edge2_freq(:));

% Enable window
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'on');
end

% Update handles structure
guidata(hObject, handles);

%% Draw Campbell diagram
function DrawCampbell(handles,undock)

% Plot settings and axis limits
if undock
    Plot = figure();
    set(Plot, 'Name', 'Campbell diagram')
else
    axis(handles.Campbell);
end

cla reset;
set(gca, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'Fontsize', 8);
xlabel('Operating range [rpm]')
ylabel('Mode frequency [Hz]')

% Axis limits
RPM1 = handles.Control.Torque.SpeedB/handles.Drivetrain.Gearbox.Ratio;
RPM2 = handles.Control.Torque.SpeedC/handles.Drivetrain.Gearbox.Ratio;
RPMlim = handles.RPMs(end);
Flim = [];
if get(handles.Check_Tower_ForeAft1,'Value')
    Flim = [Flim; max(handles.Tower_ForeAft1_freq)];
end
if get(handles.Check_Tower_ForeAft2,'Value')
    Flim = [Flim; max(handles.Tower_ForeAft2_freq)];
end
if get(handles.Check_Tower_SideSide1,'Value')
    Flim = [Flim; max(handles.Tower_SideSide1_freq)];
end
if get(handles.Check_Tower_SideSide2,'Value')
    Flim = [Flim; max(handles.Tower_SideSide2_freq)];
end
if get(handles.Check_Blade_Flap1,'Value')
    Flim = [Flim; max(handles.Blade_Flap1_freq)];
end
if get(handles.Check_Blade_Flap2,'Value')
    Flim = [Flim; max(handles.Blade_Flap2_freq)];
end
if get(handles.Check_Blade_Edge1,'Value')
    Flim = [Flim; max(handles.Blade_Edge1_freq)];
end
if get(handles.Check_Blade_Edge2,'Value')
    Flim = [Flim; max(handles.Blade_Edge2_freq)];
end
if isempty(Flim) || sum(isnan(Flim)) > 0
    Flim = 1;
end
Flim = ceil(max(Flim));
xlim([0 RPMlim])
ylim([0 Flim(1)])

% Plot 1-12P lines
hold on
patch([RPM1 RPM2 RPM2 RPM1], [0 0 1e3 1e3], [0.95 0.95 0.95], 'EdgeColor', 'none');
plot([1 1]*RPM1, [0 1e3], '--', 'Color', [0.5 0.5 0.5])
plot([1 1]*RPM2, [0 1e3], '--', 'Color', [0.5 0.5 0.5])
for P = 1:12
    eval(['draw = get(handles.Check_', int2str(P), 'P, ''Value'');']);
    if draw
        plot([0 RPMlim], P*[0 RPMlim]/60, 'Color', [0.25 0.25 0.25])
        if P*0.95*RPMlim/60 > 0.95*Flim
            x = 0.95*Flim*60/P;
            y = 0.95*Flim;
            if x > RPM1 && x < RPM2
                bg = [0.95 0.95 0.95];
            else
                bg = 'w';
            end
            text(x, y, [int2str(P), 'P'], 'Fontsize', 8, 'Color', [0.25 0.25 0.25], 'BackgroundColor', bg, 'HorizontalAlignment', 'center')
        else
            x = 0.95*RPMlim;
            y = P*0.95*RPMlim/60;
            if x > RPM1 && x < RPM2
                bg = [0.95 0.95 0.95];
            else
                bg = 'w';
            end
            text(x, y, [int2str(P), 'P'], 'Fontsize', 8, 'Color', [0.25 0.25 0.25], 'BackgroundColor', bg, 'HorizontalAlignment', 'center')
        end
    end
end

% Plot natural frequencies
LegendEntries = {};
Lines = [];
i = 1;
if get(handles.Check_Tower_ForeAft1,'Value')
    l1 = plot(handles.RPMs, handles.Tower_ForeAft1_freq, '-o', 'Color', [127.5 0 0]/255, 'MarkerSize', 2);
    LegendEntries{i} = '1st tower fore-aft';
    Lines(i) = l1;
    i = i + 1;
end
if get(handles.Check_Tower_ForeAft2,'Value')
    l2 = plot(handles.RPMs, handles.Tower_ForeAft2_freq, '-o', 'Color', [237 28 36]/255, 'MarkerSize', 2);
    LegendEntries{i} = '2nd tower fore-aft';
    Lines(i) = l2;
    i = i + 1;
end
if get(handles.Check_Tower_SideSide1,'Value')
    l3 = plot(handles.RPMs, handles.Tower_SideSide1_freq, '-o', 'Color',  [241 140 34]/255, 'MarkerSize', 2);
    LegendEntries{i} = '1st tower side-side';
    Lines(i) = l3;
    i = i + 1;
end
if get(handles.Check_Tower_SideSide2,'Value')
    l4 = plot(handles.RPMs, handles.Tower_SideSide2_freq, '-o', 'Color', [255 222 23]/255, 'MarkerSize', 2);
    LegendEntries{i} = '2nd tower side-side';
    Lines(i) = l4;
    i = i + 1;
end
if get(handles.Check_Blade_Flap1,'Value')
    l5 = plot(handles.RPMs, handles.Blade_Flap1_freq, '-o', 'Color', [173 209 54]/255, 'MarkerSize', 2);
    LegendEntries{i} = '1st blade flapwise';
    Lines(i) = l5;
    i = i + 1;
end
if get(handles.Check_Blade_Flap2,'Value')
    l6 = plot(handles.RPMs, handles.Blade_Flap2_freq, '-o', 'Color', [8 135 67]/255, 'MarkerSize', 2);
    LegendEntries{i} = '2nd blade flapwise';
    Lines(i) = l6;
    i = i + 1;
end
if get(handles.Check_Blade_Edge1,'Value')
    l7 = plot(handles.RPMs, handles.Blade_Edge1_freq, '-o', 'Color', [71 195 211]/255, 'MarkerSize', 2);
    LegendEntries{i} = '1st blade edgewise';
    Lines(i) = l7;
    i = i + 1;
end
if get(handles.Check_Blade_Edge2,'Value')
    l8 = plot(handles.RPMs, handles.Blade_Edge2_freq, '-o', 'Color',  [33 64 154]/255, 'MarkerSize', 2);
    LegendEntries{i} = '2nd blade edgewise';
    Lines(i) = l8;
end
hold off

if ~isempty(Lines)
    legend(Lines, LegendEntries, 'Location', 'NorthWest', 'Color', 'none', 'Box', 'off');
end

%% Undock button
function Undock_Callback(hObject, eventdata, handles)
DrawCampbell(handles,true)

%% Plot tower mode shapes
function PlotTowerMode(shape,mode,handles)

% Get geometry from handles
Tower = handles.Tower;

% Create fixed number of tower elements (N = 21)
N = 21;
x = linspace(Tower.Height(1), Tower.Height(end), N);
Tower.Diameter = interp1(Tower.Height, Tower.Diameter, x);
Tower.Height = x;

% Plot settings
Plot = figure();
set(Plot, 'Name', mode)
view(135,20)
light
lightangle(0, 45)
hold on
axis equal
axis off

% Undeformed tower geometry
r = repmat(Tower.Diameter(:)/2,[1,N]);
azi = repmat(linspace(0,2*pi,N),[length(Tower.Diameter),1]);
x = r.*cos(azi);
y = r.*sin(azi);
z = repmat(Tower.Height(:),[1,N]);
surf(x,y,z, ...
    'FaceColor', [1, 1, 1], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit', ...
    'FaceAlpha', 0.25)

% Wall
vert = [[-1 1 1 -1 -1 1 1 -1] * Tower.Diameter(1); ...
        [-1 -1 1 1 -1 -1 1 1] * Tower.Diameter(1); ...
        [-1 -1 -1 -1 0 0 0 0]/10];
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices', vert', 'Faces', faces, ...
    'FaceColor', [0.75, 0.75, 0.75], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit')

% Deformed tower geometry
N = 20;
r = repmat(Tower.Diameter(:)/2,[1,N]);
azi = repmat(linspace(0,2*pi,N),[length(Tower.Diameter),1]);
x = r.*cos(azi) + repmat(Tower.Height(end)*shape(1,:)',[1,N]);
y = r.*sin(azi) + repmat(Tower.Height(end)*shape(2,:)',[1,N]);
z = repmat(Tower.Height(:),[1,N]);
surf(x,y,z, ...
    'FaceColor', [1, 1, 1], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit')

%% Plot blade mode shapes
function PlotBladeMode(shape,mode,handles)

% Get geometry from handles
Blade = handles.Blade;
Airfoil = handles.Airfoil;

% Create fixed number of blade elements (N = 21)
N = 21;
x = linspace(Blade.Radius(1), Blade.Radius(end), N);
Blade.Twist = interp1(Blade.Radius, Blade.Twist, x);
Blade.Chord = interp1(Blade.Radius, Blade.Chord, x);
Blade.PitchAxis = interp1(Blade.Radius, Blade.PitchAxis, x);
Blade.Thickness = interp1(Blade.Radius, Blade.Thickness, x);
Blade.NFoil = interp1(Blade.Radius, Blade.NFoil, x, 'nearest');
Blade.Radius = x;

% Plot settings
Plot = figure();
set(Plot, 'Name', mode)
view(45,20)
light
lightangle(0, 45)
hold on
axis equal
axis off

% Undeformed blade geometry
N = 40;
x = [];
y = [];
for i = 1:length(Blade.NFoil)
    x = [x; ([handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(1,1:N:end), handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(1,1)] - Blade.PitchAxis(i)) * Blade.Chord(i)];
end
for i = 1:length(Blade.NFoil)
    t_u = max(handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,1:200));
    t_l = min(handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,200:end));
    t(i) = t_u - t_l;
end
for i = 1:length(Blade.NFoil)
    y = [y; -1*[handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,1:N:end), handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,1)] * Blade.Thickness(i)/t(i)];
end
z = repmat(Blade.Radius(:),[1,size(x,2)]);
for i = 1:length(Blade.NFoil)
    t = pi/2 - Blade.Twist(i) * pi/180;
    Rz = [cos(t),-sin(t), 0; ...
          sin(t), cos(t), 0; ...
          0,      0,      1];
    A = Rz * [x(i,:); y(i,:); z(i,:)];
    x(i,:) = A(1,:);
    y(i,:) = A(2,:);
    z(i,:) = A(3,:);
end
surf(z,y,x, ...
    'FaceColor', [1, 1, 1], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit', ...
    'FaceAlpha', 0.25)

% Wall
vert = [[-1 -1 -1 -1 0 0 0 0]/10 + Blade.Radius(1); ...
        [-1 -1 1 1 -1 -1 1 1] * Blade.Chord(1); ...
        [-1 1 1 -1 -1 1 1 -1] * Blade.Chord(1)];
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
patch('Vertices', vert', 'Faces', faces, ...
    'FaceColor', [0.75, 0.75, 0.75], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit')

% Deformed blade geometry
N = 40;
x = [];
y = [];
for i = 1:length(Blade.NFoil)
    x = [x; ([handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(1,1:N:end), handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(1,1)] - Blade.PitchAxis(i)) * Blade.Chord(i)];
end
for i = 1:length(Blade.NFoil)
    t_u = max(handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,1:200));
    t_l = min(handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,200:end));
    t(i) = t_u - t_l;
end
for i = 1:length(Blade.NFoil)
    y = [y; -1*[handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,1:N:end), handles.Airfoil.Geometry{Blade.IFoil(Blade.NFoil(i))}(2,1)] * Blade.Thickness(i)/t(i)];
end
z = repmat(Blade.Radius(:),[1,size(x,2)]);
for i = 1:length(Blade.NFoil)
    t = pi/2 - Blade.Twist(i) * 2*pi/180 - pi/180*shape(3,i);
    Rz = [cos(t),-sin(t), 0; ...
          sin(t), cos(t), 0; ...
          0,      0,      1];
    A = Rz * [x(i,:); y(i,:); z(i,:)];
    x(i,:) = A(1,:) + Blade.Radius(end)*shape(1,i);
    y(i,:) = A(2,:) + Blade.Radius(end)*shape(2,i);
    z(i,:) = A(3,:);
end
surf(z,y,x, ...
    'FaceColor', [1, 1, 1], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit')

%% Frequency output boxes (inactive)
function Tower_ForeAft1_textbox_Callback(hObject, eventdata, handles)
function Tower_ForeAft1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Tower_ForeAft2_textbox_Callback(hObject, eventdata, handles)
function Tower_ForeAft2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Tower_SideSide1_textbox_Callback(hObject, eventdata, handles)
function Tower_SideSide1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Tower_SideSide2_textbox_Callback(hObject, eventdata, handles)
function Tower_SideSide2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Blade_Flap1_textbox_Callback(hObject, eventdata, handles)
function Blade_Flap1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Blade_Flap2_textbox_Callback(hObject, eventdata, handles)
function Blade_Flap2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Blade_Edge1_textbox_Callback(hObject, eventdata, handles)
function Blade_Edge1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Blade_Edge2_textbox_Callback(hObject, eventdata, handles)
function Blade_Edge2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Plot buttons
function Plot_Tower_ForeAft1_Callback(hObject, eventdata, handles)
PlotTowerMode(handles.Tower_ForeAft1_shape, ['First tower fore-aft mode (', num2str(handles.Tower_ForeAft1_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Tower_ForeAft2_Callback(hObject, eventdata, handles)
PlotTowerMode(handles.Tower_ForeAft2_shape, ['Second tower fore-aft mode (', num2str(handles.Tower_ForeAft2_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Tower_SideSide1_Callback(hObject, eventdata, handles)
PlotTowerMode(handles.Tower_SideSide1_shape, ['First tower side-side mode (', num2str(handles.Tower_SideSide1_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Tower_SideSide2_Callback(hObject, eventdata, handles)
PlotTowerMode(handles.Tower_SideSide2_shape, ['Second tower side-side mode (', num2str(handles.Tower_SideSide2_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Blade_Flap1_Callback(hObject, eventdata, handles)
PlotBladeMode(handles.Blade_Flap1_shape, ['First blade flapwise mode (', num2str(handles.Blade_Flap1_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Blade_Flap2_Callback(hObject, eventdata, handles)
PlotBladeMode(handles.Blade_Flap2_shape, ['Second blade flapwise mode (', num2str(handles.Blade_Flap2_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Blade_Edge1_Callback(hObject, eventdata, handles)
PlotBladeMode(handles.Blade_Edge1_shape, ['First blade edgewise mode (', num2str(handles.Blade_Edge1_freq, '%0.2f'), ' Hz)'],handles)
function Plot_Blade_Edge2_Callback(hObject, eventdata, handles)
PlotBladeMode(handles.Blade_Edge2_shape, ['Second blade edgewise mode (', num2str(handles.Blade_Edge2_freq, '%0.2f'), ' Hz)'],handles)

%% Rotation speed step size - textbox
function RPMstep_textbox_Callback(hObject, eventdata, handles)
RPMend = handles.Control.Torque.SpeedC/handles.Drivetrain.Gearbox.Ratio;
RPMend = ceil(RPMend/5)*5;
RPMstep = str2double(get(handles.RPMstep_textbox, 'String'));
handles.RPMs = 0:RPMstep:RPMend;
guidata(hObject, handles);
function RPMstep_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% NP frequencies - check boxes
function Check_1P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_2P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_3P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_4P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_5P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_6P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_7P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_8P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_9P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_10P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_11P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_12P_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)

%% Mode frequencies - check boxes
function Check_Blade_Edge2_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Blade_Edge1_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Blade_Flap2_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Blade_Flap1_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function checkbox20_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Tower_SideSide2_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Tower_SideSide1_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Tower_ForeAft2_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
function Check_Tower_ForeAft1_Callback(hObject, eventdata, handles)
DrawCampbell(handles,false)
