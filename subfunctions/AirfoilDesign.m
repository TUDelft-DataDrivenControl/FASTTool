%% Initialization code
function varargout = AirfoilDesign(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AirfoilDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @AirfoilDesign_OutputFcn, ...
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
function AirfoilDesign_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Airfoil = varargin{1};
handles.Blade = varargin{2};
handles.Input = varargin;

% Update input fields
set(handles.AirfoilCollection, 'String', handles.Airfoil.Name)

% Selected airfoil
i = get(handles.AirfoilCollection,'Value');
handles.CurrentAirfoil = i;

% Update stall characteristics
set(handles.StallAngle1_textbox,'String',num2str(handles.Airfoil.StallAngle1(i)));
set(handles.StallAngle2_textbox,'String',num2str(handles.Airfoil.StallAngle2(i)));
set(handles.CnSlope_textbox,'String',num2str(handles.Airfoil.CnSlope(i)));
set(handles.CritCn1_textbox,'String',num2str(handles.Airfoil.CritCn1(i)));
set(handles.CritCn2_textbox,'String',num2str(handles.Airfoil.CritCn2(i)));

% Update tables
set(handles.AirfoilPerformance_Table, 'Data', ...
   [num2cell(handles.Airfoil.Alpha{i}), ...
    num2cell(handles.Airfoil.Cl{i}), ...
    num2cell(handles.Airfoil.Cd{i}), ...
    num2cell(handles.Airfoil.Cm{i})]);
set(handles.AirfoilCoordinates_Table, 'Data', ...
   [num2cell(handles.Airfoil.Geometry{i}(1,:)'), ...
    num2cell(handles.Airfoil.Geometry{i}(2,:)')]);
set(handles.AirfoilPerformance_TableSize_textbox, 'String', length(handles.Airfoil.Alpha{i}));
set(handles.AirfoilPerformance_TableSize_slider, 'Value', length(handles.Airfoil.Alpha{i}));
set(handles.AirfoilCoordinates_TableSize_textbox, 'String', length(handles.Airfoil.Geometry{i}(1,:)));
set(handles.AirfoilCoordinates_TableSize_slider, 'Value', length(handles.Airfoil.Geometry{i}(1,:)));

% Update airfoil plot
DrawAirfoil(handles)

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.AirfoilDesign);

%% Closing function
function AirfoilDesign_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.AirfoilDesign);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.AirfoilDesign);

%% Output function
function varargout = AirfoilDesign_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save
    
    % Get table data
    AirfoilPerformance_Table = get(handles.AirfoilPerformance_Table, 'Data');
    AirfoilCoordinates_Table = get(handles.AirfoilCoordinates_Table, 'Data');

    % Find invalid cells
    for i = 1:size(AirfoilPerformance_Table,1)
        for j = 1:size(AirfoilPerformance_Table,2)
            invalid(i,j) = ...
                isempty(cell2mat(AirfoilPerformance_Table(i,j))) + ...
                sum(isnan(cell2mat(AirfoilPerformance_Table(i,j))));
        end
    end

    % Fill vectors
    Alpha = [];
    Cl = [];
    Cd = [];
    Cm = [];
    for i = 1:size(AirfoilPerformance_Table,1)
        if sum(invalid(i,:)) == 0
            Alpha = [Alpha; AirfoilPerformance_Table{i,1}];
            Cl = [Cl; AirfoilPerformance_Table{i,2}];
            Cd = [Cd; AirfoilPerformance_Table{i,3}];
            Cm = [Cm; AirfoilPerformance_Table{i,4}];
        end
    end

    % Find invalid cells
    for i = 1:size(AirfoilCoordinates_Table,1)
        for j = 1:size(AirfoilCoordinates_Table,2)
            invalid(i,j) = ...
                isempty(cell2mat(AirfoilCoordinates_Table(i,j))) + ...
                sum(isnan(cell2mat(AirfoilCoordinates_Table(i,j))));
        end
    end

    % Fill vectors
    Geometry = [];
    for i = 1:size(AirfoilCoordinates_Table,1)
        if sum(invalid(i,:)) == 0
            Geometry = [Geometry, cell2mat(AirfoilCoordinates_Table(i,:))'];
        end
    end
    
    % Update values
    i = handles.CurrentAirfoil;
    handles.Airfoil.Alpha{i} = Alpha;
    handles.Airfoil.Cl{i} = Cl;
    handles.Airfoil.Cd{i} = Cd;
    handles.Airfoil.Cm{i} = Cm;
    handles.Airfoil.Geometry{i} = Geometry;
    handles.Airfoil.StallAngle1(i) = str2double(get(handles.StallAngle1_textbox,'String'));
    handles.Airfoil.StallAngle2(i) = str2double(get(handles.StallAngle2_textbox,'String'));
    handles.Airfoil.CnSlope(i) = str2double(get(handles.CnSlope_textbox,'String'));
    handles.Airfoil.CritCn1(i) = str2double(get(handles.CritCn1_textbox,'String'));
    handles.Airfoil.CritCn2(i) = str2double(get(handles.CritCn2_textbox,'String'));
    
    % Refine geometry
    for i = 1:length(handles.Airfoil.Geometry)
        
        % Airfoil coordinates
        x = handles.Airfoil.Geometry{i}(1,:);
        y = handles.Airfoil.Geometry{i}(2,:);
        
        % Find upper and lower surface
        j = find(x == min(x));
        j = j(round(length(j)/2));
        if mean(y(1:j)) > 0
            x_u = x(1:j);
            x_l = x(j:end);
            y_u = y(1:j);
            y_l = y(j:end);
        else
            x_u = x(j:end);
            x_l = x(1:j);
            y_u = y(j:end);
            y_l = y(1:j);
        end
        
        % Remove duplicates
        [x_u,j] = unique(x_u);
        y_u = y_u(j);
        [x_l,j] = unique(x_l);
        y_l = y_l(j);
        
        % Normalize chord length
        x_u = x_u - min(x_u);
        x_u = x_u/max(x_u);
        x_l = x_l - min(x_l);
        x_l = x_l/max(x_l);
        
        % Interpolate
        x = 0.5 + 0.5*cos(linspace(0,2*pi,399));
        y_u = interp1(x_u, y_u, x(1:200), 'pchip');
        y_l = interp1(x_l, y_l, x(200:end), 'pchip');
        y = [y_u(1:200), y_l(2:end)];
        
        % Update geometry
        handles.Airfoil.Geometry{i} = zeros(2,399);
        handles.Airfoil.Geometry{i}(1,:) = x;
        handles.Airfoil.Geometry{i}(2,:) = y;
        
        % Complete airfoil data with flat plate behavior of the NACA 64-618
        if i > 2
            if min(handles.Airfoil.Alpha{i}) > -180
                handles.Airfoil.Cl{i} = [handles.Airfoil.Cl{8}(handles.Airfoil.Alpha{8} < min(handles.Airfoil.Alpha{i})); handles.Airfoil.Cl{i}];
                handles.Airfoil.Cd{i} = [handles.Airfoil.Cd{8}(handles.Airfoil.Alpha{8} < min(handles.Airfoil.Alpha{i})); handles.Airfoil.Cd{i}];
                handles.Airfoil.Cm{i} = [handles.Airfoil.Cm{8}(handles.Airfoil.Alpha{8} < min(handles.Airfoil.Alpha{i})); handles.Airfoil.Cm{i}];
                handles.Airfoil.Alpha{i} = [handles.Airfoil.Alpha{8}(handles.Airfoil.Alpha{8} < min(handles.Airfoil.Alpha{i})); handles.Airfoil.Alpha{i}];
            end
            if max(handles.Airfoil.Alpha{i}) < 180
                handles.Airfoil.Cl{i} = [handles.Airfoil.Cl{i}; handles.Airfoil.Cl{8}(handles.Airfoil.Alpha{8} > max(handles.Airfoil.Alpha{i}))];
                handles.Airfoil.Cd{i} = [handles.Airfoil.Cd{i}; handles.Airfoil.Cd{8}(handles.Airfoil.Alpha{8} > max(handles.Airfoil.Alpha{i}))];
                handles.Airfoil.Cm{i} = [handles.Airfoil.Cm{i}; handles.Airfoil.Cm{8}(handles.Airfoil.Alpha{8} > max(handles.Airfoil.Alpha{i}))];
                handles.Airfoil.Alpha{i} = [handles.Airfoil.Alpha{i}; handles.Airfoil.Alpha{8}(handles.Airfoil.Alpha{8} > max(handles.Airfoil.Alpha{i}))];
            end
        end
    end
    
    varargout{1} = handles.Airfoil;
    varargout{2} = handles.Blade;
    
else
    
    varargout = handles.Input;
    
end

% Close figure
delete(hObject)

%% New airfoil
function NewAirfoil_Callback(hObject, eventdata, handles)

% Prompt for new airfoil name
name = inputdlg('Enter airfoil name:');

if ~isempty(name)

    % Make room in airfoil structure
    i = length(handles.Airfoil.Name) + 1;
    handles.Airfoil.Name{i} = name{1};
    handles.Airfoil.Geometry{i} = [1,0,1;0,0,0];
    handles.Airfoil.Alpha{i} = [-180;0;180];
    handles.Airfoil.Cl{i} = [0;0;0];
    handles.Airfoil.Cd{i} = [0;0;0];
    handles.Airfoil.Cm{i} = [0;0;0];
    handles.Airfoil.StallAngle1(i) = 0;
    handles.Airfoil.StallAngle2(i) = 0;
    handles.Airfoil.CnSlope(i) = 0;
    handles.Airfoil.CritCn1(i) = 0;
    handles.Airfoil.CritCn2(i) = 0;

    % Update list
    list = get(handles.AirfoilCollection,'String');
    list{i} = name{1};
    set(handles.AirfoilCollection,'String',list);
    set(handles.AirfoilCollection,'Value',i);

    % Update handles structure
    guidata(hObject, handles);
    AirfoilCollection_Callback(handles.AirfoilCollection, eventdata, handles);

    % Update airfoil plot
    DrawAirfoil(handles)

end

%% Delete airfoil
function DeleteAirfoil_Callback(hObject, eventdata, handles)

% Prevent deletion of base airfoils
if length(handles.Airfoil.Name) <= 8
    errordlg('Cannot delete base airfoils.');
else
    
    % Delete airfoil
    i = get(handles.AirfoilCollection,'Value');
    handles.Airfoil.Name(i) = [];
    handles.Airfoil.Geometry(i) = [];
    handles.Airfoil.Alpha(i) = [];
    handles.Airfoil.Cl(i) = [];
    handles.Airfoil.Cd(i) = [];
    handles.Airfoil.Cm(i) = [];
    handles.Airfoil.StallAngle1(i) = [];
    handles.Airfoil.StallAngle2(i) = [];
    handles.Airfoil.CnSlope(i) = [];
    handles.Airfoil.CritCn1(i) = [];
    handles.Airfoil.CritCn2(i) = [];
    handles.Blade.IFoil(handles.Blade.IFoil == i) = i - 1;
    
    % Update list
    set(handles.AirfoilCollection,'String',handles.Airfoil.Name);
    if i > length(handles.Airfoil.Name)
        set(handles.AirfoilCollection,'Value',i-1);
    end
    
    % Update handles structure
    guidata(hObject, handles);
    AirfoilCollection_Callback(handles.AirfoilCollection, eventdata, handles);
        
end

%% Airfoil selection
function AirfoilCollection_Callback(hObject, eventdata, handles)

% Get table data
AirfoilPerformance_Table = get(handles.AirfoilPerformance_Table, 'Data');
AirfoilCoordinates_Table = get(handles.AirfoilCoordinates_Table, 'Data');

% Find invalid cells
for i = 1:size(AirfoilPerformance_Table,1)
    for j = 1:size(AirfoilPerformance_Table,2)
        invalid(i,j) = ...
            isempty(cell2mat(AirfoilPerformance_Table(i,j))) + ...
            sum(isnan(cell2mat(AirfoilPerformance_Table(i,j))));
    end
end

% Fill vectors
Alpha = [];
Cl = [];
Cd = [];
Cm = [];
for i = 1:size(AirfoilPerformance_Table,1)
    if sum(invalid(i,:)) == 0
        Alpha = [Alpha; AirfoilPerformance_Table{i,1}];
        Cl = [Cl; AirfoilPerformance_Table{i,2}];
        Cd = [Cd; AirfoilPerformance_Table{i,3}];
        Cm = [Cm; AirfoilPerformance_Table{i,4}];
    end
end

% Find invalid cells
for i = 1:size(AirfoilCoordinates_Table,1)
    for j = 1:size(AirfoilCoordinates_Table,2)
        invalid(i,j) = ...
            isempty(cell2mat(AirfoilCoordinates_Table(i,j))) + ...
            sum(isnan(cell2mat(AirfoilCoordinates_Table(i,j))));
    end
end

% Fill vectors
Geometry = [];
for i = 1:size(AirfoilCoordinates_Table,1)
    if sum(invalid(i,:)) == 0
        Geometry = [Geometry, cell2mat(AirfoilCoordinates_Table(i,:))'];
    end
end

% Update values
i = handles.CurrentAirfoil;
handles.Airfoil.Alpha{i} = Alpha;
handles.Airfoil.Cl{i} = Cl;
handles.Airfoil.Cd{i} = Cd;
handles.Airfoil.Cm{i} = Cm;
handles.Airfoil.Geometry{i} = Geometry;
handles.Airfoil.StallAngle1(i) = str2double(get(handles.StallAngle1_textbox,'String'));
handles.Airfoil.StallAngle2(i) = str2double(get(handles.StallAngle2_textbox,'String'));
handles.Airfoil.CnSlope(i) = str2double(get(handles.CnSlope_textbox,'String'));
handles.Airfoil.CritCn1(i) = str2double(get(handles.CritCn1_textbox,'String'));
handles.Airfoil.CritCn2(i) = str2double(get(handles.CritCn2_textbox,'String'));

% Selected airfoil
i = get(hObject,'Value');
handles.CurrentAirfoil = i;

% Update stall characteristics
set(handles.StallAngle1_textbox,'String',num2str(handles.Airfoil.StallAngle1(i)));
set(handles.StallAngle2_textbox,'String',num2str(handles.Airfoil.StallAngle2(i)));
set(handles.CnSlope_textbox,'String',num2str(handles.Airfoil.CnSlope(i)));
set(handles.CritCn1_textbox,'String',num2str(handles.Airfoil.CritCn1(i)));
set(handles.CritCn2_textbox,'String',num2str(handles.Airfoil.CritCn2(i)));

% Update tables
set(handles.AirfoilPerformance_Table, 'Data', ...
   [num2cell(handles.Airfoil.Alpha{i}), ...
    num2cell(handles.Airfoil.Cl{i}), ...
    num2cell(handles.Airfoil.Cd{i}), ...
    num2cell(handles.Airfoil.Cm{i})]);
set(handles.AirfoilCoordinates_Table, 'Data', ...
   [num2cell(handles.Airfoil.Geometry{i}(1,:)'), ...
    num2cell(handles.Airfoil.Geometry{i}(2,:)')]);
set(handles.AirfoilPerformance_TableSize_textbox, 'String', length(handles.Airfoil.Alpha{i}));
set(handles.AirfoilPerformance_TableSize_slider, 'Value', length(handles.Airfoil.Alpha{i}));
set(handles.AirfoilCoordinates_TableSize_textbox, 'String', length(handles.Airfoil.Geometry{i}(1,:)));
set(handles.AirfoilCoordinates_TableSize_slider, 'Value', length(handles.Airfoil.Geometry{i}(1,:)));

% Update airfoil plot
DrawAirfoil(handles)

% Update handles structure
guidata(hObject, handles);
function AirfoilCollection_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Plot airfoil polars
function PlotButton_Callback(hObject, eventdata, handles)

% Plot settings
Plot = figure();
set(Plot, 'Name', 'Airfoil polars')

% Lift curve
subplot(2,1,1)
hold on
plot([-180 180], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.25)
alpha = cell2mat(handles.AirfoilPerformance_Table.Data(:,1));
Cl = cell2mat(handles.AirfoilPerformance_Table.Data(:,2));
Cd = cell2mat(handles.AirfoilPerformance_Table.Data(:,3));
Cm = cell2mat(handles.AirfoilPerformance_Table.Data(:,4));
Cl_plot = plot(alpha, Cl);
Cd_plot = plot(alpha, Cd);
Cm_plot = plot(alpha, Cm);
set(gca, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'XTick', -180:30:180, ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'Fontsize', 8);
xlim([-180 180])
legend([Cl_plot, Cd_plot, Cm_plot], {'Cl'; 'Cd'; 'Cm'}, 'Location', 'Best', 'Color', 'none', 'Box', 'off');
xlabel('Angle of attack [deg]')
ylabel('[-]')

subplot(2,1,2)
hold on
Cl = Cl(alpha > -30);
Cd = Cd(alpha > -30);
alpha = alpha(alpha > -30);
Cl = Cl(alpha < 30);
Cd = Cd(alpha < 30);
alpha = alpha(alpha < 30);

hold on
plot([0 0.025], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.25)
plot(Cd,Cl)
hold off
set(gca, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'Fontsize', 8);
xlim([0 0.025])
xlabel('Cd [-]')
ylabel('Cl [-]')

%% Draw airfoil
function DrawAirfoil(handles)

% Plot settings
axis(handles.AirfoilPlot);
cla reset;
set(gca, ...
    'XTick', 0.0:0.1:1.0, ...
    'YTick', -0.5:0.1:0.5, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Box', 'on', ...
    'Layer', 'top', ...
    'Fontsize', 8);
axis equal
xlim([-0.1 1.1])
ylim([-0.6 0.6])

% Plot shape
x = cell2mat(handles.AirfoilCoordinates_Table.Data(:,1));
y = cell2mat(handles.AirfoilCoordinates_Table.Data(:,2));
hold on
plot([-0.1 1.1], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.25)
patch(x, y, [0 165, 213]/255, 'EdgeColor', 0.5*[0 165, 213]/255, 'LineWidth', 1)
hold off

%% Airfoil performance - table size text box
function AirfoilPerformance_TableSize_textbox_Callback(hObject, eventdata, handles)

rows = ceil(str2double(get(hObject, 'String')));
Table = get(handles.AirfoilPerformance_Table, 'Data');
if rows < 2
    rows = 2;
elseif rows > get(handles.AirfoilPerformance_TableSize_slider, 'Max')
    rows = get(handles.AirfoilPerformance_TableSize_slider, 'Max');
elseif isnan(rows)
    rows = size(Table,1);
end
if rows < size(Table,1)
    Table = Table(1:rows,:);
    set(handles.AirfoilPerformance_Table, 'Data', Table);
elseif rows > size(Table,1)
    Table{rows,size(Table,2)} = [];
    set(handles.AirfoilPerformance_Table, 'Data', Table);
end
set(hObject, 'String', int2str(rows));
set(handles.AirfoilPerformance_TableSize_slider, 'Value', rows);
function AirfoilPerformance_TableSize_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Airfoil performance - table size slider
function AirfoilPerformance_TableSize_slider_Callback(hObject, eventdata, handles)
rows = get(hObject, 'Value');
Table = get(handles.AirfoilPerformance_Table, 'Data');
if rows < 2
    set(hObject, 'Value', 2);
else
    if rows < size(Table,1)
        Table = Table(1:rows,:);
        set(handles.AirfoilPerformance_Table, 'Data', Table);
    elseif rows > size(Table,1)
        Table{rows,size(Table,2)} = [];
        set(handles.AirfoilPerformance_Table, 'Data', Table);
    end
    set(handles.AirfoilPerformance_TableSize_textbox, 'String', int2str(rows));
end
function AirfoilPerformance_TableSize_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Airfoil performance - edit cells checkbox
function AirfoilPerformance_EditCells_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1
    set(handles.AirfoilPerformance_Table, 'ColumnEditable', [true true true true]);
else
    set(handles.AirfoilPerformance_Table, 'ColumnEditable', [false false false false]);
end

%% Airfoil performance - cell selection
function AirfoilPerformance_Table_CellSelectionCallback(hObject, eventdata, handles)
handles.Selection = eventdata.Indices;
guidata(hObject, handles); 

%% Airfoil performance - keyboard functions
function AirfoilPerformance_Table_KeyPressFcn(hObject, eventdata, handles)

% Paste functionality
if strcmpi(char(eventdata.Modifier),'control') && strcmp(eventdata.Key, 'v')
    
    % Get and reshape clipboard data
    Paste = clipboard('paste');
    rows = numel(strsplit(strtrim(Paste), '\n'));
    Paste = strsplit(strtrim(Paste));
    columns = numel(Paste)/rows;
    Paste = reshape(Paste,columns,rows)';
    
    % Convert numeric cells
    for i = 1:rows
        for j = 1:columns
        	Paste{i,j} = str2double(Paste{i,j});
        end
    end

    % Target cells
    Table = get(hObject, 'Data');
    Selection = handles.Selection;
    i1 = Selection(1,1);
    i2 = Selection(1,1) + (rows-1);
    j1 = Selection(1,2);
    j2 = Selection(1,2) + (columns-1);
    if i2 > size(Table,1)
        if i2 > get(handles.AirfoilPerformance_TableSize_slider, 'Max')
            i2 = get(handles.AirfoilPerformance_TableSize_slider, 'Max');
        end
        set(handles.AirfoilPerformance_TableSize_textbox, 'String', int2str(i2));
        set(handles.AirfoilPerformance_TableSize_slider, 'Value', i2);
    end
    if j2 > size(Table,2)
        j2 = size(Table,2);
    end
    
    % Constrain data types
    Table(i1:i2,j1:j2) = Paste(1:(1+i2-i1),1:(1+j2-j1));
    for i = i1:i2
        for j = 1:size(Table,2)
            if ~isnumeric(Table{i,j})
                Table{i,j} = NaN;
            end
        end
    end

    % Update airfoilperformance_table
    set(hObject, 'Data', Table);

end

%% Airfoil coordinates - table size text box
function AirfoilCoordinates_TableSize_textbox_Callback(hObject, eventdata, handles)

rows = ceil(str2double(get(hObject, 'String')));
Table = get(handles.AirfoilCoordinates_Table, 'Data');
if rows < 2
    rows = 2;
elseif rows > get(handles.AirfoilCoordinates_TableSize_slider, 'Max')
    rows = get(handles.AirfoilCoordinates_TableSize_slider, 'Max');
elseif isnan(rows)
    rows = size(Table,1);
end
if rows < size(Table,1)
    Table = Table(1:rows,:);
    set(handles.AirfoilCoordinates_Table, 'Data', Table);
elseif rows > size(Table,1)
    Table{rows,size(Table,2)} = [];
    set(handles.AirfoilCoordinates_Table, 'Data', Table);
end
set(hObject, 'String', int2str(rows));
set(handles.AirfoilCoordinates_TableSize_slider, 'Value', rows);
function AirfoilCoordinates_TableSize_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Airfoil coordinates - table size slider
function AirfoilCoordinates_TableSize_slider_Callback(hObject, eventdata, handles)
rows = get(hObject, 'Value');
Table = get(handles.AirfoilCoordinates_Table, 'Data');
if rows < 2
    set(hObject, 'Value', 2);
else
    if rows < size(Table,1)
        Table = Table(1:rows,:);
        set(handles.AirfoilCoordinates_Table, 'Data', Table);
    elseif rows > size(Table,1)
        Table{rows,size(Table,2)} = [];
        set(handles.AirfoilCoordinates_Table, 'Data', Table);
    end
    set(handles.AirfoilCoordinates_TableSize_textbox, 'String', int2str(rows));
end
function AirfoilCoordinates_TableSize_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Airfoil coordinates - edit cells checkbox
function AirfoilCoordinates_EditCells_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1
    set(handles.AirfoilCoordinates_Table, 'ColumnEditable', [true true]);
else
    set(handles.AirfoilCoordinates_Table, 'ColumnEditable', [false false]);
end

%% Airfoil coordinates - cell selection
function AirfoilCoordinates_Table_CellSelectionCallback(hObject, eventdata, handles)
handles.Selection = eventdata.Indices;
guidata(hObject, handles);

%% Airfoil coordinates - keyboard functions
function AirfoilCoordinates_Table_KeyPressFcn(hObject, eventdata, handles)

% Paste functionality
if strcmpi(char(eventdata.Modifier),'control') && strcmp(eventdata.Key, 'v')
    
    % Get and reshape clipboard data
    Paste = clipboard('paste');
    rows = numel(strsplit(strtrim(Paste), '\n'));
    Paste = strsplit(strtrim(Paste));
    columns = numel(Paste)/rows;
    Paste = reshape(Paste,columns,rows)';
    
    % Convert numeric cells
    for i = 1:rows
        for j = 1:columns
        	Paste{i,j} = str2double(Paste{i,j});
        end
    end

    % Target cells
    Table = get(hObject, 'Data');
    Selection = handles.Selection;
    i1 = Selection(1,1);
    i2 = Selection(1,1) + (rows-1);
    j1 = Selection(1,2);
    j2 = Selection(1,2) + (columns-1);
    if i2 > size(Table,1)
        if i2 > get(handles.AirfoilCoordinates_TableSize_slider, 'Max')
            i2 = get(handles.AirfoilCoordinates_TableSize_slider, 'Max');
        end
        set(handles.AirfoilCoordinates_TableSize_textbox, 'String', int2str(i2));
        set(handles.AirfoilCoordinates_TableSize_slider, 'Value', i2);
    end
    if j2 > size(Table,2)
        j2 = size(Table,2);
    end
    
    % Constrain data types
    Table(i1:i2,j1:j2) = Paste(1:(1+i2-i1),1:(1+j2-j1));
    for i = i1:i2
        for j = 1:size(Table,2)
            if ~isnumeric(Table{i,j})
                Table{i,j} = NaN;
            end
        end
    end

    % Update airfoilperformance_table
    set(hObject, 'Data', Table);

end

%% Airfoil coordinates - cell edit redraw
function AirfoilCoordinates_Table_CellEditCallback(hObject, eventdata, handles)

% Selected airfoil
i = get(handles.AirfoilCollection,'Value');

% Update airfoil plot
DrawAirfoil(handles)

%% Positive stall angle - text box
function StallAngle1_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.StallAngle1))
end
handles.StallAngle1 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function StallAngle1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Negative stall angle - text box
function StallAngle2_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.StallAngle2))
end
handles.StallAngle2 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function StallAngle2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Cn slope - text box
function CnSlope_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CnSlope))
end
handles.CnSlope = str2double(get(hObject,'String'));
guidata(hObject, handles);
function CnSlope_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Critical Cn at positive AoA - text box
function CritCn1_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CritCn1))
end
handles.CritCn1 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function CritCn1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Critical Cn at negative AoA - text box
function CritCn2_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CritCn2))
end
handles.CritCn2 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function CritCn2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
