%% Initialization code
function varargout = TowerDesign(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TowerDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @TowerDesign_OutputFcn, ...
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
function TowerDesign_OpeningFcn(hObject, eventdata, handles, varargin)

% Set background image
h = axes('Units', 'Normalized', 'position', [0 0 1 1]);
uistack(h, 'bottom');
img = imread('graphics\tower.png');
imagesc(img);
set(h, 'HandleVisibility', 'off', 'visible','off')

% Get input
handles.Tower = varargin{1};
handles.TowerStyle = varargin{2};
handles.TowerPaint = varargin{3};
handles.Input = varargin(1);

% Update input fields
set(handles.TableSize_textbox, 'String', length(handles.Tower.Height));
set(handles.TableSize_slider, 'Value', length(handles.Tower.Height));
set(handles.Tower_BottomThickness_textbox, 'String', num2str(handles.Tower.BottomThickness));
set(handles.Tower_TopThickness_textbox, 'String', num2str(handles.Tower.TopThickness));
set(handles.Hub_Height_textbox, 'String', num2str(handles.Tower.HubHeight));
set(handles.Table, 'Data', ...
   [num2cell(handles.Tower.Height(:)), ...
    num2cell(handles.Tower.Diameter(:)), ...
    num2cell(handles.Tower.ExtraMass(:))]);

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.TowerDesign);

%% Closing function
function TowerDesign_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.TowerDesign);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.TowerDesign);

%% Output function
function varargout = TowerDesign_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save

    % Get geometry from table
    Table = get(handles.Table, 'Data');

    % Find invalid cells
    for i = 1:size(Table,1)
        for j = 1:size(Table,2)
            invalid(i,j) = ...
                isempty(cell2mat(Table(i,j))) + ...
                sum(isnan(cell2mat(Table(i,j))));
        end
    end

    % Extract geometry from table
    Tower.Height = [];
    Tower.Diameter = [];
    Tower.ExtraMass = [];
    for i = 1:size(Table,1)
        if sum(invalid(i,:)) == 0
            Tower.Height = [Tower.Height; Table{i,1}];
            Tower.Diameter = [Tower.Diameter; Table{i,2}];
            Tower.ExtraMass = [Tower.ExtraMass; Table{i,3}];
        end
    end
    Tower.BottomThickness = str2double(get(handles.Tower_BottomThickness_textbox, 'String'));
    Tower.TopThickness = str2double(get(handles.Tower_TopThickness_textbox, 'String'));
    Tower.HubHeight = str2double(get(handles.Hub_Height_textbox, 'String'));
    
    % Derived properties
    Tower.ShearModulus = 80.8e9;
    Tower.YoungsModulus = 210e9;
    Tower.Density = 7850;
    Tower.WallThickness = linspace(Tower.BottomThickness,Tower.TopThickness,length(Tower.Height))';
    Tower.Mass = Tower.Density * pi* (Tower.Diameter.^2 - (Tower.Diameter-2*Tower.WallThickness).^2) ./ 4 + Tower.ExtraMass;
    Tower.EI = Tower.YoungsModulus * pi/64*(Tower.Diameter.^4 - (Tower.Diameter-2*Tower.WallThickness).^4);
    Tower.GJ = Tower.ShearModulus * pi/32*(Tower.Diameter.^4 - (Tower.Diameter-2*Tower.WallThickness).^4);
    Tower.EA = Tower.YoungsModulus * pi*Tower.Diameter.*Tower.WallThickness;
    Tower.Iner = 0.5 * Tower.Mass .* (Tower.Diameter/2 - Tower.WallThickness).^2;
    
    varargout{1} = Tower;
    
else
    
    varargout = handles.Input;
    
end

% Close figure
delete(hObject)

%% Hub height - text box
function Hub_Height_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Tower.HubHeight))
end
handles.Tower.HubHeight = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Hub_Height_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Tower top wall thickness - text box
function Tower_TopThickness_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Tower.TopThickness))
end
handles.Tower.TopThickness = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Tower_TopThickness_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Tower bottom wall thickness - text box
function Tower_BottomThickness_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Tower.BottomThickness))
end
handles.Tower.BottomThickness = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Tower_BottomThickness_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Table size - text box
function TableSize_textbox_Callback(hObject, eventdata, handles)

rows = ceil(str2double(get(hObject, 'String')));
Table = get(handles.Table, 'Data');
if rows < 2
    rows = 2;
elseif rows > get(handles.TableSize_slider, 'Max')
    rows = get(handles.TableSize_slider, 'Max');
elseif isnan(rows)
    rows = size(Table,1);
end
if rows < size(Table,1)
    Table = Table(1:rows,:);
    set(handles.Table, 'Data', Table);
elseif rows > size(Table,1)
    Table{rows,size(Table,2)} = [];
    set(handles.Table, 'Data', Table);
end
set(hObject, 'String', int2str(rows));
set(handles.TableSize_slider, 'Value', rows);
function TableSize_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Table size - slider
function TableSize_slider_Callback(hObject, eventdata, handles)
rows = get(hObject, 'Value');
Table = get(handles.Table, 'Data');
if rows < 2
    set(hObject, 'Value', 2);
else
    if rows < size(Table,1)
        Table = Table(1:rows,:);
        set(handles.Table, 'Data', Table);
    elseif rows > size(Table,1)
        Table{rows,size(Table,2)} = [];
        set(handles.Table, 'Data', Table);
    end
    set(handles.TableSize_textbox, 'String', int2str(rows));
end
function TableSize_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Edit cells - checkbox
function EditCells_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1
    set(handles.Table, 'ColumnEditable', [true true true true true true true]);
else
    set(handles.Table, 'ColumnEditable', [false false false false false false false]);
end

%% Table cell selection
function Table_CellSelectionCallback(hObject, eventdata, handles)
handles.Selection = eventdata.Indices;
guidata(hObject, handles); 

%% Table keyboard functions
function Table_KeyPressFcn(hObject, eventdata, handles)

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
        if i2 > get(handles.TableSize_slider, 'Max')
            i2 = get(handles.TableSize_slider, 'Max');
        end
        set(handles.TableSize_textbox, 'String', int2str(i2));
        set(handles.TableSize_slider, 'Value', i2);
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

    % Update table
    set(hObject, 'Data', Table);

end

%% Plot geometry
function PlotButton_Callback(hObject, eventdata, handles)

% Extract geometry from table
Table = get(handles.Table, 'Data');
for i = 1:size(Table,1)
    if isempty(Table{i,1})
        Tower.Height(i) = NaN;
    else
        Tower.Height(i) = Table{i,1};
    end
    if isempty(Table{i,2})
        Tower.Diameter(i) = NaN;
    else
        Tower.Diameter(i) = Table{i,2};
    end
end
Tower.BottomThickness = str2double(get(handles.Tower_BottomThickness_textbox, 'String'));
Tower.TopThickness = str2double(get(handles.Tower_TopThickness_textbox, 'String'));
Tower_WallThickness = linspace(Tower.BottomThickness,Tower.TopThickness,length(Tower.Height));

% Color scheme
EdgeColor = 'none';
White = [240, 240, 240]/255;
Grey = [120, 120, 120]/255;

% Set axis
Plot = figure();
set(Plot, 'Name', 'Tower plot')
view(0,180/pi*atan(sin(pi/4)))
light
lightangle(0, 45)
axis equal
axis off
xlim([-0.5 0.5] * max(Tower.Diameter));
ylim([-0.5 0.5] * max(Tower.Diameter));
zlim([min(Tower.Height) max(Tower.Height)]);
hold on

% Plot tower
N = 200 + 1;
s = Tower.Height/Tower.Height(end);
if handles.TowerStyle == 1
    paint = ones(size(s));
elseif handles.TowerStyle == 2
    paint = zeros(size(s));
    s = 0.05*round(20*s);
    paint(s >= 0.2) = 1;
    paint(s >= 0.4) = 0;
    paint(s >= 0.6) = 1;
    paint(s >= 0.8) = 0;
    paint(s >= 1.0) = 1;
elseif handles.TowerStyle == 3
    paint = zeros(size(s));
    paint(s > 0.25) = 1;
    paint(s > 0.31) = 0;
elseif handles.TowerStyle == 4
    paint = (s > 0.8);
elseif handles.TowerStyle == 5
    paint = zeros(size(s));
    paint(Tower.Height <= 30) = 1;
elseif handles.TowerStyle == 6
    paint = zeros(size(s));
    n = length(s(Tower.Height <= 40));
    paint(1:n) = linspace(1,0,n);
elseif handles.TowerStyle == 7
    paint = [ones(1,size(s,2)); zeros(1,size(s,2))];
    s = 0.05*round(20*s);
    paint(1,s >= 0.1) = 0;
    paint(1,s >= 0.2) = 1;
    paint(1,s >= 0.3) = 0;
    paint(1,s >= 0.4) = 1;
    paint(1,s >= 0.5) = 0;
    paint(1,s >= 0.6) = 1;
    paint(1,s >= 0.7) = 0;
    paint(1,s >= 0.8) = 1;
    paint(1,s >= 0.9) = 0;
    paint(1,s >= 1.0) = 1;
    paint(2,s >= 0.1) = 1;
    paint(2,s >= 0.2) = 0;
    paint(2,s >= 0.3) = 1;
    paint(2,s >= 0.4) = 0;
    paint(2,s >= 0.5) = 1;
    paint(2,s >= 0.6) = 0;
    paint(2,s >= 0.7) = 1;
    paint(2,s >= 0.8) = 0;
    paint(2,s >= 0.9) = 1;
    paint(2,s >= 1.0) = 0;
    paint = [paint; paint];
elseif handles.TowerStyle
    paint = (s < 0.15);
end
if handles.TowerStyle == 7
    paint = kron(paint,ones(round(N/4),1));
    paint = paint(1:N-1,:);
    TowerColors = zeros([size(paint),3]);
    TowerColors(:,:,1) = paint*handles.TowerPaint(1) + (1-paint)*White(1);
    TowerColors(:,:,2) = paint*handles.TowerPaint(2) + (1-paint)*White(2);
    TowerColors(:,:,3) = paint*handles.TowerPaint(3) + (1-paint)*White(3);
    TowerColors = permute(TowerColors,[2,1,3]);
else
    TowerColors = kron(paint(:),handles.TowerPaint) + kron((1-paint(:)),White);
    TowerColors = repmat(TowerColors,[1,1,N]);
    TowerColors = permute(TowerColors,[1,3,2]);
end

r = repmat(Tower.Diameter(:)/2,[1,N]);
azi = repmat(linspace(0,2*pi,N),[length(Tower.Diameter),1]);
x = r.*cos(azi);
y = r.*sin(azi);
z = repmat(Tower.Height(:),[1,N]);

TowerColors = [TowerColors, TowerColors(:,end,:)];
surf(x,y,z, ...
    'CData', TowerColors, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', EdgeColor, ...
    'AmbientStrength', 0.5, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.5, ...
    'BackFaceLighting', 'reverselit')

for i = 1:length(Tower.Height)
    x_(i,:) = x(i,:) * (Tower.Diameter(i) - 2*Tower_WallThickness(i))/Tower.Diameter(i);
    y_(i,:) = y(i,:) * (Tower.Diameter(i) - 2*Tower_WallThickness(i))/Tower.Diameter(i);
end
z_ = z;

surf(x_,y_,z_, ...
    'FaceColor', Grey, ...
    'EdgeColor', EdgeColor, ...
    'AmbientStrength', 0.5, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.5, ...
    'BackFaceLighting', 'reverselit')

patch(...
    [x(1,:), flip(x_(1,:))], ...
    [y(1,:), flip(y_(1,:))], ...
    [z(1,:), flip(z_(1,:))], ...
    'g', 'FaceColor', White, ...
    'EdgeColor', EdgeColor, ...
    'AmbientStrength', 0.5, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.5, ...
    'BackFaceLighting', 'reverselit')

patch(...
    [x(end,:), flip(x_(end,:))], ...
    [y(end,:), flip(y_(end,:))], ...
    [z(end,:), flip(z_(end,:))], ...
    'g', 'FaceColor', White, ...
    'EdgeColor', EdgeColor, ...
    'AmbientStrength', 0.5, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.5, ...
    'BackFaceLighting', 'reverselit')
