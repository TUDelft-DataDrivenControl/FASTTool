%% Initialization code
function varargout = PitchGain(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PitchGain_OpeningFcn, ...
                   'gui_OutputFcn',  @PitchGain_OutputFcn, ...
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
function PitchGain_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Control = varargin{1};
handles.Input = varargin;

% Update input fields
set(handles.LPFCutOff_text, 'String', num2str(handles.Control.Pitch.LowPassCutOffFreq));

% TODO: get rid of the field checking and update the baseline's MAT 
% file in the future release
if ~isfield(handles.Control.Pitch, 'Notch2_beta1GS') % executed once
    handles.Control.Pitch.Notch2_beta1GS = zeros(14,1);
end
if ~isfield(handles.Control.Pitch, 'Notch2_beta2GS') % executed once
    handles.Control.Pitch.Notch2_beta2GS = zeros(14,1);
end
if ~isfield(handles.Control.Pitch, 'Notch2_wnGS') % executed once
    handles.Control.Pitch.Notch2_wnGS = zeros(14,1);
end

handles.TableSize = length(handles.Control.Pitch.KpGS);
for i = 1:handles.TableSize
    if isnan(handles.Control.Pitch.ScheduledPitchAngles(i))
        TableData(i,1) = num2cell([]);
    else
        TableData(i,1) = num2cell(handles.Control.Pitch.ScheduledPitchAngles(i)*180/pi);
    end
    if isnan(handles.Control.Pitch.KpGS(i))
        TableData(i,2) = num2cell([]);
    else
        TableData(i,2) = num2cell(handles.Control.Pitch.KpGS(i));
    end
    if isnan(handles.Control.Pitch.KiGS(i))
        TableData(i,3) = num2cell([]);
    else
        TableData(i,3) = num2cell(handles.Control.Pitch.KiGS(i));
    end
    if isnan(handles.Control.Pitch.Notch_beta1GS(i))
        TableData(i,4) = num2cell([]);
    else
        TableData(i,4) = num2cell(handles.Control.Pitch.Notch_beta1GS(i));
    end
    if isnan(handles.Control.Pitch.Notch_beta2GS(i))
        TableData(i,5) = num2cell([]);
    else
        TableData(i,5) = num2cell(handles.Control.Pitch.Notch_beta2GS(i));
    end
    if isnan(handles.Control.Pitch.Notch_wnGS(i))
        TableData(i,6) = num2cell([]);
    else
        TableData(i,6) = num2cell(handles.Control.Pitch.Notch_wnGS(i));
    end
    % --- Second notch filter start
    if isnan(handles.Control.Pitch.Notch2_beta1GS(i))
        TableData(i,7) = num2cell([]);
    else
        TableData(i,7) = num2cell(handles.Control.Pitch.Notch2_beta1GS(i));
    end
    
    if isnan(handles.Control.Pitch.Notch2_beta2GS(i))
        TableData(i,8) = num2cell([]);
    else
        TableData(i,8) = num2cell(handles.Control.Pitch.Notch2_beta2GS(i));
    end
    
    if isnan(handles.Control.Pitch.Notch2_wnGS(i))
        TableData(i,9) = num2cell([]);
    else
        TableData(i,9) = num2cell(handles.Control.Pitch.Notch2_wnGS(i));
    end
    % --- Second notch filter end
    if isnan(handles.Control.Pitch.LowPassCutOffFreqGS(i))
        TableData(i,10) = num2cell([]);
    else
        TableData(i,10) = num2cell(handles.Control.Pitch.LowPassCutOffFreqGS(i));
    end
end
set(handles.Constant_Ki_textbox, 'String', num2str(handles.Control.Pitch.Ki));
set(handles.Constant_Kp_textbox, 'String', num2str(handles.Control.Pitch.Kp));

set(handles.Constant_Notch_B1_textbox, 'String', num2str(handles.Control.Pitch.Notch_beta1));
set(handles.Constant_Notch_B2_textbox, 'String', num2str(handles.Control.Pitch.Notch_beta2));
set(handles.Constant_Notch_wn_textbox, 'String', num2str(handles.Control.Pitch.Notch_wn));
% --- Second notch filter start
% TODO: get rid of the field checking and update the baseline's MAT 
% file in the future release
if ~isfield(handles.Control.Pitch, 'Notch2_beta1')
    handles.Control.Pitch.Notch2_beta1 = 0;
end
if ~isfield(handles.Control.Pitch, 'Notch2_beta2')
    handles.Control.Pitch.Notch2_beta2 = 0;
end
if ~isfield(handles.Control.Pitch, 'Notch2_wn')
    handles.Control.Pitch.Notch2_wn = 0;
end
set(handles.Constant_Notch2_B1_textbox, 'String', num2str(handles.Control.Pitch.Notch2_beta1));
set(handles.Constant_Notch2_B2_textbox, 'String', num2str(handles.Control.Pitch.Notch2_beta2));
set(handles.Constant_Notch2_wn_textbox, 'String', num2str(handles.Control.Pitch.Notch2_wn));
% --- Second notch filter end

set(handles.ExportPlotData_pushbutton, 'Enable', 'off')
set(handles.RefreshPlot_pushbutton, 'Enable', 'off')


set(handles.TableSize_textbox, 'String', length(handles.Control.Pitch.KpGS));
set(handles.TableSize_slider, 'Value', length(handles.Control.Pitch.KpGS));
set(handles.Table, 'Data', TableData);

% Check for gain scheduling
if handles.Control.Pitch.Scheduled
    set(handles.ConstantGain, 'Value', 0);
    set(handles.GainScheduled, 'Value', 1);
    set(handles.Constant_Ki_textbox, 'Enable', 'off')
    set(handles.Constant_Kp_textbox, 'Enable', 'off')
    set(handles.LPFCutOff_text, 'Enable', 'off')
    set(handles.Constant_Notch_B1_textbox, 'Enable', 'off')
    set(handles.Constant_Notch_B2_textbox, 'Enable', 'off')
    set(handles.Constant_Notch_wn_textbox, 'Enable', 'off')
    set(handles.Constant_Notch2_B1_textbox, 'Enable', 'off')
    set(handles.Constant_Notch2_B2_textbox, 'Enable', 'off')
    set(handles.Constant_Notch2_wn_textbox, 'Enable', 'off')
    set(handles.Table, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Value', 0)
    set(handles.TableSize_textbox, 'Enable', 'off')
    set(handles.TableSize_slider, 'Enable', 'off')
else
    set(handles.ConstantGain, 'Value', 1);
    set(handles.GainScheduled, 'Value', 0);
    set(handles.Constant_Ki_textbox, 'Enable', 'on')
    set(handles.Constant_Kp_textbox, 'Enable', 'on')
    set(handles.LPFCutOff_text, 'Enable', 'on')
    set(handles.Constant_Notch_B1_textbox, 'Enable', 'on')
    set(handles.Constant_Notch_B2_textbox, 'Enable', 'on')
    set(handles.Constant_Notch_wn_textbox, 'Enable', 'on')
    set(handles.Constant_Notch2_B1_textbox, 'Enable', 'on')
    set(handles.Constant_Notch2_B2_textbox, 'Enable', 'on')
    set(handles.Constant_Notch2_wn_textbox, 'Enable', 'on')
    set(handles.Table, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Value', 0)
    set(handles.TableSize_textbox, 'Enable', 'off')
    set(handles.TableSize_slider, 'Enable', 'off')
end

if handles.Control.Pitch.LowPassOrder == 1
    set(handles.FirstOrderLPF_radio, 'Value', 1);
    set(handles.SecondOrderLPF_radio, 'Value', 0);
else
    set(handles.FirstOrderLPF_radio, 'Value', 0);
    set(handles.SecondOrderLPF_radio, 'Value', 1);
end

% Update handles structure
guidata(hObject, handles);

% Adapt to screen size
monitor = get(groot, 'MonitorPosition');
width = monitor(1,3);
window = get(handles.Window_group, 'Position');
needed_width = window(3)+16;
slider_pos = get(handles.Window_slider, 'Position');
pos = get(hObject, 'Position');
if width >= needed_width
    set(hObject, 'Position', [pos(1) pos(2) pos(3) pos(4)-slider_pos(4)]);
    set(handles.Window_group, 'Position', [0 0 window(3) window(4)]);
    set(handles.Window_slider, 'Visible', 'off');
else
    set(handles.Window_slider, 'Position', [0 0 pos(3) slider_pos(4)]);
end

% Halt window
uiwait(handles.PitchGain);

%% Closing function
function PitchGain_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.PitchGain);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.PitchGain);

%% Output function
function varargout = PitchGain_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save

    handles = UpdateHandlesWithTableData(handles);
    
    handles.Control.Pitch.Ki = str2double(get(handles.Constant_Ki_textbox,'String'));
    handles.Control.Pitch.Kp = str2double(get(handles.Constant_Kp_textbox,'String'));
    
    handles.Control.Pitch.Notch_beta1 = str2double(get(handles.Constant_Notch_B1_textbox,'String'));
    handles.Control.Pitch.Notch_beta2 = str2double(get(handles.Constant_Notch_B2_textbox,'String'));
    handles.Control.Pitch.Notch_wn = str2double(get(handles.Constant_Notch_wn_textbox,'String'));

    % Second notch start
    handles.Control.Pitch.Notch2_beta1 = str2double(get(handles.Constant_Notch2_B1_textbox,'String'));
    handles.Control.Pitch.Notch2_beta2 = str2double(get(handles.Constant_Notch2_B2_textbox,'String'));
    handles.Control.Pitch.Notch2_wn = str2double(get(handles.Constant_Notch2_wn_textbox,'String'));
    % Second notch end

    if get(handles.FirstOrderLPF_radio, 'Value')
        handles.Control.Pitch.LowPassOrder = 1;
    else
        handles.Control.Pitch.LowPassOrder = 2;
    end
    
    handles.Control.Pitch.Scheduled = get(handles.GainScheduled,'Value');
    
    % Do a check if the controller parameters make sense
    assertControllerGains(handles)
    
    varargout{1} = handles.Control;
else
    varargout = handles.Input;
end

% Close figure
delete(hObject)

%% Edit cells - checkbox
function EditCells_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1
    set(handles.Table, 'ColumnEditable', true(1, handles.TableSize));
    set(handles.TableSize_textbox, 'Enable', 'on')
    set(handles.TableSize_slider, 'Enable', 'on')
else
    set(handles.Table, 'ColumnEditable', false(1, handles.TableSize));
    set(handles.TableSize_textbox, 'Enable', 'off')
    set(handles.TableSize_slider, 'Enable', 'off')
end

%% Table cell selection
function Table_CellSelectionCallback(hObject, eventdata, handles)
handles.Selection = eventdata.Indices;
guidata(hObject, handles); 

%% Table keyboard functions
function Table_KeyPressFcn(hObject, eventdata, handles)

% Paste functionality
if strcmpi(char(eventdata.Modifier),'control') && strcmp(eventdata.Key, 'v') && all(get(handles.Table, 'ColumnEditable'))
    
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
        if i2 > handles.TableSize
            i2 = handles.TableSize;
        end
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

%% Constant gain - radio button
function ConstantGain_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.Constant_Ki_textbox, 'Enable', 'on')
    set(handles.Constant_Kp_textbox, 'Enable', 'on')
    set(handles.Constant_Notch_B1_textbox, 'Enable', 'on')
    set(handles.Constant_Notch_B2_textbox, 'Enable', 'on')
    set(handles.Constant_Notch_wn_textbox, 'Enable', 'on')
    set(handles.Constant_Notch2_B1_textbox, 'Enable', 'on')
    set(handles.Constant_Notch2_B2_textbox, 'Enable', 'on')
    set(handles.Constant_Notch2_wn_textbox, 'Enable', 'on')
    set(handles.LPFCutOff_text, 'Enable', 'on')
    set(handles.Table, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Enable', 'off')
    
    handles.Control.Pitch.Scheduled = get(handles.GainScheduled,'Value');
    guidata(hObject, handles);
end
    
%% Gain scheduling - radio button
function GainScheduled_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.Constant_Ki_textbox, 'Enable', 'off')
    set(handles.Constant_Kp_textbox, 'Enable', 'off')
    set(handles.Constant_Notch_B1_textbox, 'Enable', 'off')
    set(handles.Constant_Notch_B2_textbox, 'Enable', 'off')
    set(handles.Constant_Notch_wn_textbox, 'Enable', 'off')
    set(handles.Constant_Notch2_B1_textbox, 'Enable', 'off')
    set(handles.Constant_Notch2_B2_textbox, 'Enable', 'off')
    set(handles.Constant_Notch2_wn_textbox, 'Enable', 'off')
    set(handles.LPFCutOff_text, 'Enable', 'off')
    set(handles.Table, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Enable', 'on')
   
    handles.Control.Pitch.Scheduled = get(handles.GainScheduled,'Value');
    guidata(hObject, handles);
end

%% Constant integration constant - text box
function Constant_Ki_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Ki(end)))
end
handles.Control.Pitch.Ki(end) = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Ki_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Constant gain - text box
function Constant_Kp_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Kp(end)))
end
handles.Control.Pitch.Kp(end) = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Kp_textbox_CreateFcn(hObject, eventdata, handles)
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

%% --- Executes on button press in checkboxes.
function PlotLPF_checkbox_Callback(hObject, eventdata, handles)
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    if AllDisabled
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableButtons(handles, 'on')
    end

    BodePlot(handles, false, false)

function PlotPI_checkbox_Callback(hObject, eventdata, handles)
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    if AllDisabled
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableButtons(handles, 'on')
    end
	
    BodePlot(handles, false, false)


function PlotNotch_checkbox_Callback(hObject, eventdata, handles)
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    if AllDisabled
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableButtons(handles, 'on')
    end

    BodePlot(handles, false, false)

function PlotNom_checkbox_Callback(hObject, eventdata, handles)
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    if AllDisabled
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableButtons(handles, 'on')
    end

    BodePlot(handles, false, false)

function PlotLoopGain_checkbox_Callback(hObject, eventdata, handles)
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    if AllDisabled
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableButtons(handles, 'on')
    end
    
	BodePlot(handles, false, false)

function PlotGMPM_checkbox_Callback(hObject, eventdata, handles)
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    if AllDisabled
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableButtons(handles, 'on')
    end
    
	BodePlot(handles, false, false)

%% --- Executes on button presses
function UndockBode_pushbutton_Callback(hObject, eventdata, handles)
    BodePlot(handles, true, false)

function BodePlot(handles, undock, exportData)
   	cla(handles.BodeMag_axes,'reset')
    cla(handles.BodePhase_axes,'reset')
    clc; % to clear all warnings
    
    % Create vector of plot colors
    plotCol = linspace(0.8, 0, length(handles.SelectedListboxContents));
    plotLineStyle = {'-', '--', '-'};
    plotLineWidth = [1, 1, 2];

    % Evaluate state of buttons checked
    [AllDisabled, ~] = CheckStateCheckboxes(handles);
    
    % Create transfer function of filters in series according to selection GUI
    Plant = tf(1,1)*ones(1,length(handles.SelectedListboxContents));
    for i = 1:length(handles.SelectedListboxContents)
        selIndex = findnearest(str2double(handles.SelectedListboxContents{i}), handles.Lin.Pitch*180/pi);
        iGenSpeed = find(contains(handles.sysm{selIndex,1}.OutputName, 'ED GenSpeed'));
        iBlPitchCPC = find(contains(handles.sysm{selIndex,1}.InputName, 'collective blade-pitch'));
        Plant(1,i) = pi/30*handles.sysm{selIndex,1}(iGenSpeed,iBlPitchCPC);
    end
    
    % Synthesize controller for plotting the loopgain
    ControllerLG = calculateController(handles, true);
    
    % Synthesize controller for plotting the controller only
    Controller = calculateController(handles, false);

    % Evaluate the user choise to only see the controller, or to show the
    % loop-gain
    LoopGain = tf(1,1)*ones(1,length(handles.SelectedListboxContents));
    for i = 1:length(handles.SelectedListboxContents)
        LoopGain(1,i) = series(ControllerLG(:,i), Plant(:,i));
    end
    
    w_llimit = -2; 
    w_ulimit = 2;
    w = logspace(w_llimit, w_ulimit, 1000);

    [~, Controller_MagResponse, Controller_PhaseResponse] = ...
        calculateFreqResp(Controller, w);
    
    [LoopGainFRF, LoopGainMagResponse, LoopGainPhaseResponse] = ...
        calculateFreqResp(LoopGain, w);
    
    [PlantFRF, PlantMagResponse, PlantPhaseResponse] = ...
        calculateFreqResp(Plant, w);

    % Ensure the data is always saved column-wise for single line plots
    if size(Controller_MagResponse, 1) == 1
        Controller_MagResponse = Controller_MagResponse(:);
        Controller_PhaseResponse = Controller_PhaseResponse(:);
        LoopGainMagResponse = LoopGainMagResponse(:);
        LoopGainPhaseResponse = LoopGainPhaseResponse(:);
        PlantMagResponse = PlantMagResponse(:);
        PlantPhaseResponse = PlantPhaseResponse(:);
    end
    
    % Stability margins calculation
    cell_prealloc = cell(1,length(handles.SelectedListboxContents));
    S = struct('GainMargin', cell_prealloc, ...
        'GMFrequency', cell_prealloc, ...
        'PhaseMargin', cell_prealloc, ...
        'PMFrequency', cell_prealloc, ...
        'DelayMargin', cell_prealloc, ...
        'DMFrequency', cell_prealloc, ...
        'Stable', cell_prealloc);
    GM = nan(1, length(handles.SelectedListboxContents)); 
    PM = GM;
    GMFreq = GM;
    PMFreq = GM;
    IsGMFreqWithinLimit = true(1, length(handles.SelectedListboxContents));
    IsPMFreqWithinLimit = IsGMFreqWithinLimit;
    warning('off','backtrace') % suppress warning
    if get(handles.PlotGMPM_checkbox,'value')
        LoopGainMagResponseAbs = db2mag(LoopGainMagResponse); % for margins calculation
        for i = 1:length(handles.SelectedListboxContents)    
            S(i) = allmargin(LoopGainMagResponseAbs(:,i),LoopGainPhaseResponse(:,i),w');
            try
                GM(i) = mag2db(S(i).GainMargin(1));
                PM(i) = S(i).PhaseMargin(1);
            catch
                S(i) = allmargin(LoopGain(1,i));
                try
                    GM(i) = mag2db(S(i).GainMargin(1));
                    PM(i) = S(i).PhaseMargin(1);
                catch
                    quest_msg = sprintf("Cannot obtain GM and/or PM for the model with pitch angle %.3f [deg]! Please try to increase the gain (Kp or Ki) and try again!", ...
                        string(handles.SelectedListboxContents{i}));
                    questdlg(quest_msg, 'Error', 'OK', 'OK');
                    continue;
                end
            end
            
            msg = "";
            GMFreq(i) = S(i).GMFrequency(1);
            IsGMFreqWithinLimit(i) = (log10(GMFreq(i)) >= w_llimit) && (log10(GMFreq(i)) <= w_ulimit);
            if ~IsGMFreqWithinLimit(i)
                msg = sprintf("Gain margin frequency is outside the figure axis limit: %.5f! [rad/s] ", ...
                    GMFreq(i));
            end
            
            PMFreq(i) = S(i).PMFrequency(1);
            IsPMFreqWithinLimit(i) = (log10(PMFreq(i)) >= w_llimit) && (log10(PMFreq(i)) <= w_ulimit);
            if ~IsPMFreqWithinLimit(i)
                msg = msg + sprintf("Phase margin frequency is outside the figure axis limit: %.5f! [rad/s] ", ...
                    PMFreq(i));
            end

            if ~IsGMFreqWithinLimit(i) || ~IsPMFreqWithinLimit(i)
                msg = "Cannot draw one or more stability margins! " + msg;
                warning(msg);
            end
        end
    end
    warning('on','backtrace')
    
    if get(handles.PlotGMPM_checkbox,'value')
        legendMargins = cell(length(handles.SelectedListboxContents), 12);
        legendMargins(:,1) = {' GM: '};
        legendMargins(:,3) = {' [dB] '};
        legendMargins(:,4) = {' \omega_{pc}: '};
        legendMargins(:,6) = {' [rad/s] '};
        legendMargins(:,7) = {' PM: '};
        legendMargins(:,9) = {' [deg] '};
        legendMargins(:,10) = {' \omega_{gc}: '};
        legendMargins(:,12) = {' [rad/s] '}; 
        
        for i = 1:length(handles.SelectedListboxContents)    
            legendMargins{i,2} = num2str(GM(i),'%.2f');
            legendMargins{i,5} = num2str(GMFreq(i),'%.4f');
            legendMargins{i,8} = num2str(PM(i), '%.2f');
            legendMargins{i,11} = num2str(PMFreq(i),'%.4f');
        end
        legendMargins = strcat(legendMargins(:,1), legendMargins(:,2), ...
            legendMargins(:,3), legendMargins(:,4), legendMargins(:,5), ...
            legendMargins(:,6), legendMargins(:,7), legendMargins(:,8), ...
            legendMargins(:,9), legendMargins(:,10), legendMargins(:,11), ...
            legendMargins(:,12));
    end
    
    if exportData
        ControllerLG_FRF = freqresp(ControllerLG, w);
        
        frd_Plant = frd(PlantFRF, w);
        frd_Controller = frd(ControllerLG_FRF, w);
        frd_LoopGain = frd(LoopGainFRF, w);
        
        PitchAngles = str2double(handles.SelectedListboxContents);
        
        variablesToSave = {'frd_Plant', 'frd_Controller', 'frd_LoopGain', 'PitchAngles'};
        
        if get(handles.PlotGMPM_checkbox,'value')
            variablesToSave = [variablesToSave {'GM', 'PM', 'GMFreq', 'PMFreq'}];
        end
        
        uisave(variablesToSave, 'ExportPlotData');
    end

    if undock
        Plot = figure();
        set(Plot, 'Name', 'Bode diagram')
        subplot(211)
    else
        axes(handles.BodeMag_axes);
    end
    
    for i = 1:length(handles.SelectedListboxContents)
        if get(handles.PlotNom_checkbox, 'Value')
            h(i) = semilogx(w, PlantMagResponse(:,i), 'Color', ones(1,3)*plotCol(i), 'LineStyle', plotLineStyle{1}, 'LineWidth', plotLineWidth(1));
            hold on
        end
        if not(AllDisabled)
            h(i) = semilogx(w, Controller_MagResponse(:,i), 'Color', ones(1,3)*plotCol(i), 'LineStyle', plotLineStyle{2}, 'LineWidth', plotLineWidth(2));
            hold on
        end
        if get(handles.PlotLoopGain_checkbox, 'Value')
            if get(handles.PlotGMPM_checkbox,'value') && IsGMFreqWithinLimit(i)
                h(i) = semilogx([GMFreq(i) GMFreq(i)], [0 -GM(i)], 'Color', [0 0.8 0], 'LineStyle', plotLineStyle{2}, 'LineWidth', plotLineWidth(1));
                h(i) = semilogx(GMFreq(i), -GM(i), 'o', 'Color', [0 0.8 0], 'LineStyle', plotLineStyle{1}, 'LineWidth', plotLineWidth(1));

            end
            h(i) = semilogx(w, LoopGainMagResponse(:,i), 'Color', ones(1,3)*plotCol(i), 'LineStyle', plotLineStyle{3}, 'LineWidth', plotLineWidth(3));
            hold on
        end
    end
    semilogx(w, zeros(1,length(w)), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5)
    xticklabels([])
    ylabel('Magnitude [dB]')
    set(gca, 'XScale', 'log')
    grid on
    
    legendEntry = handles.SelectedListboxContents;
    if get(handles.PlotGMPM_checkbox,'value') && exist('legendMargins','var')
        legendEntry = strcat(legendEntry,legendMargins);
    end
    if exist('h', 'var')
        legend(h, legendEntry, 'Location', 'SouthWest')
    end
    
    if undock
        subplot(212)
    else
        axes(handles.BodePhase_axes);
    end
    for i = 1:length(handles.SelectedListboxContents)
        if get(handles.PlotNom_checkbox, 'Value')
            semilogx(w, PlantPhaseResponse(:,i), 'Color', ones(1,3)*plotCol(i), 'LineStyle', plotLineStyle{1}, 'LineWidth', plotLineWidth(1)); 
            hold on
        end
        if not(AllDisabled)
            semilogx(w, Controller_PhaseResponse(:,i), 'Color', ones(1,3)*plotCol(i), 'LineStyle', plotLineStyle{2}, 'LineWidth', plotLineWidth(2)); 
            hold on
        end
        if get(handles.PlotLoopGain_checkbox, 'Value')
            semilogx(w, LoopGainPhaseResponse(:,i), 'Color', ones(1,3)*plotCol(i), 'LineStyle', plotLineStyle{3}, 'LineWidth', plotLineWidth(3)); 
            hold on
            if get(handles.PlotGMPM_checkbox,'value') && IsPMFreqWithinLimit(i)
                if PM(i) >= 0
                    dotPM = -180 + PM(i);
                else
                    dotPM = 180 + PM(i);
                end
                linePM = [-180 dotPM];
                h(i) = semilogx([PMFreq(i) PMFreq(i)], linePM, 'Color', [0 0.8 0], 'LineStyle', plotLineStyle{2}, 'LineWidth', plotLineWidth(1));
                h(i) = semilogx(PMFreq(i), dotPM, 'o', 'Color', [0 0.8 0], 'LineStyle', plotLineStyle{1}, 'LineWidth', plotLineWidth(1));
            end
        end
    end
    semilogx(w, 180*ones(1,length(w)), w, -180*ones(1,length(w)), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    xlabel('Frequency [rad/s]')
    ylabel('Phase [deg]')
    set(gca, 'XScale', 'log')
    grid on
    
    if undock
        figureHandle = findobj(Plot,'Type','axes','Visible','on');
        try
            linkaxes(figureHandle,'x')
        catch err_msg
            disp(err_msg.message)
        end
    end

function [frf, MagResponse, PhaseResponse] = calculateFreqResp(transfer_function, w)
    frf = freqresp(transfer_function, w);
    MagResponse = mag2db(squeeze(abs(frf)))';
    PhaseResponse = angle(squeeze(frf))'*180/pi;
    
% --- Executes on button press in LoadLinMat_pushbutton.
function LoadLinMat_pushbutton_Callback(hObject, eventdata, handles)
    % Get file name
    [FileName,PathName] = uigetfile('*.mat', 'Select existing linearized model');
    handles.LinModelFileName = [PathName, FileName];

    % Check if it is a valid model
    if FileName
        if exist(handles.LinModelFileName, 'file') == 2

            contents = whos('-file', handles.LinModelFileName);
            if ismember('Lin', {contents.name}) && ismember('sysm', {contents.name})
                % Load the linear models .mat-file
                load(handles.LinModelFileName)
                handles.Lin = Lin;
                handles.sysm = sysm;
                set(handles.LinFileCheck_text, 'String', 'Valid model loaded')
                set(handles.LinFileCheck_text, 'ForegroundColor', [0.0 1.0 0.0])

                % Show the above-rated linear model pitch angles in listbox
                if length(Lin.Pitch) == 1
                    LinListBoxItemsIndex = 1;
                else
                    LinListBoxItemsIndex = find([diff(Lin.Pitch(1:2)) (diff(Lin.Pitch) > 0)]);
                end
                handles.FirstAboveRatedLinModelIndex = LinListBoxItemsIndex(1);
                set(handles.LinWindSpeed_listbox, 'string', {Lin.Pitch(LinListBoxItemsIndex)*180/pi})
                
                if ~isfield(handles.Control.Pitch, 'Notch2_beta1GS')
                    handles.Control.Pitch.Notch2_beta1GS = zeros(14,1);
                end
                if ~isfield(handles.Control.Pitch, 'Notch2_beta2GS')
                    handles.Control.Pitch.Notch2_beta2GS = zeros(14,1);
                end
                if ~isfield(handles.Control.Pitch, 'Notch2_wnGS')
                    handles.Control.Pitch.Notch2_wnGS = zeros(14,1);
                end

                
                % TODO: Fix the table updating feature or just revert to
                % the previous method (manually updatin the table)
                
%                 if handles.TableSize ~= length(LinListBoxItemsIndex)
%                     handles.TableSize = length(LinListBoxItemsIndex);
%                 end
%                 
%                 % Update the table with the pitch angles data
%                 for i = 1:handles.TableSize
%                     if isnan(Lin.Pitch(LinListBoxItemsIndex(i)))
%                         TableData(i,1) = num2cell([]);
%                     else
%                         TableData(i,1) = num2cell(Lin.Pitch(LinListBoxItemsIndex(i))*180/pi);
%                     end
%                     if isnan(handles.Control.Pitch.KpGS(i))
%                         TableData(i,2) = num2cell([]);
%                     else
%                         TableData(i,2) = num2cell(handles.Control.Pitch.KpGS(i));
%                     end
%                     if isnan(handles.Control.Pitch.KiGS(i))
%                         TableData(i,3) = num2cell([]);
%                     else
%                         TableData(i,3) = num2cell(handles.Control.Pitch.KiGS(i));
%                     end
%                     if isnan(handles.Control.Pitch.Notch_beta1GS(i))
%                         TableData(i,4) = num2cell([]);
%                     else
%                         TableData(i,4) = num2cell(handles.Control.Pitch.Notch_beta1GS(i));
%                     end
%                     if isnan(handles.Control.Pitch.Notch_beta2GS(i))
%                         TableData(i,5) = num2cell([]);
%                     else
%                         TableData(i,5) = num2cell(handles.Control.Pitch.Notch_beta2GS(i));
%                     end
%                     if isnan(handles.Control.Pitch.Notch_wnGS(i))
%                         TableData(i,6) = num2cell([]);
%                     else
%                         TableData(i,6) = num2cell(handles.Control.Pitch.Notch_wnGS(i));
%                     end
%                     % --- Second notch filter
%                     if isnan(handles.Control.Pitch.Notch2_beta1GS(i))
%                         TableData(i,7) = num2cell([]);
%                     else
%                         TableData(i,7) = num2cell(handles.Control.Pitch.Notch2_beta1GS(i));
%                     end
%                     if isnan(handles.Control.Pitch.Notch2_beta2GS(i))
%                         TableData(i,8) = num2cell([]);
%                     else
%                         TableData(i,8) = num2cell(handles.Control.Pitch.Notch2_beta2GS(i));
%                     end
%                     if isnan(handles.Control.Pitch.Notch2_wnGS(i))
%                         TableData(i,9) = num2cell([]);
%                     else
%                         TableData(i,9) = num2cell(handles.Control.Pitch.Notch2_wnGS(i));
%                     end
%                     % ---
%                     if isnan(handles.Control.Pitch.LowPassCutOffFreqGS(i))
%                         TableData(i,10) = num2cell([]);
%                     else
%                         TableData(i,10) = num2cell(handles.Control.Pitch.LowPassCutOffFreqGS(i));
%                     end
%                 end
                
%                 set(handles.Table, 'Data', TableData);
                
                % Allow multiple pitch angle selection in listbox
                set(handles.LinWindSpeed_listbox, 'Max', 20, 'Min', 0);

                % Enable disabled GUI elements 
                EnableDisableListbox(handles, 'on')
                
                % Deselect all in listbox
                set(handles.LinWindSpeed_listbox, 'Value', []);
            else
                set(handles.LinFileCheck_text, 'String', 'Invalid model selected')
                set(handles.LinFileCheck_text, 'ForegroundColor', [1.0 0.0 0.0])
                
                % Disable GUI elements 
                EnableDisableListbox(handles, 'off')
                EnableDisableCheckBoxes(handles, 'off')
                EnableDisableButtons(handles, 'off')
            end
        end
    else
        % Say when an invalid file is selected
        set(handles.LinFileCheck_text, 'String', 'No file loaded')
        
        % Disable GUI elements 
        EnableDisableListbox(handles, 'off')
        EnableDisableCheckBoxes(handles, 'off')
        EnableDisableButtons(handles, 'off')
    end
    % Store in handles
    guidata(hObject, handles);

% --- Executes on selection change in LinWindSpeed_listbox.
function LinWindSpeed_listbox_Callback(hObject, eventdata, handles)
    handles.SelectedListboxContents = cellstr(get(hObject,'String'));
    handles.SelectedListboxIndex = get(hObject,'Value');
    handles.SelectedListboxContents = handles.SelectedListboxContents(handles.SelectedListboxIndex);
    if isempty(handles.SelectedListboxContents)
        EnableDisableCheckBoxes(handles, 'off')
        EnableDisableButtons(handles, 'off')
    else
        EnableDisableCheckBoxes(handles, 'on')
    end
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LinWindSpeed_listbox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function [index] = findnearest(array, value)
    [~,index] = (min(abs(array - value)));
    
function EnableDisableListbox(handles, state)
    set(handles.LinWindSpeed_listbox, 'Enable', state)
    
function EnableDisableCheckBoxes(handles, state)
    set(handles.PlotLPF_checkbox, 'Enable', state)
    set(handles.PlotPI_checkbox, 'Enable', state)
    set(handles.PlotNotch_checkbox, 'Enable', state)
    set(handles.PlotNom_checkbox, 'Enable', state)
    set(handles.PlotLoopGain_checkbox, 'Enable', state)
    set(handles.PlotGMPM_checkbox, 'Enable', state)
    
function EnableDisableButtons(handles, state)
    set(handles.UndockBode_pushbutton, 'Enable', state)
    set(handles.ExportPlotData_pushbutton, 'Enable', state)
    set(handles.RefreshPlot_pushbutton, 'Enable', state)
    
function [AllDisabled, AllControllersEnabled] = CheckStateCheckboxes(handles)
    checkBox(1) = get(handles.PlotLPF_checkbox, 'Value');
    checkBox(2) = get(handles.PlotPI_checkbox, 'Value');
    checkBox(3) = get(handles.PlotNotch_checkbox, 'Value');
    checkBox(4) = get(handles.PlotNom_checkbox, 'Value');
    checkBox(5) = get(handles.PlotLoopGain_checkbox, 'Value');
    checkBox(6) = get(handles.PlotGMPM_checkbox, 'Value');
    
    AllDisabled = all(checkBox(1:6) == 0);
    AllControllersEnabled = all(checkBox(1:3) == 1);
    
function Controller = calculateController(handles, LoopGainCheckbox)
    Controller = tf(1,1)*ones(1,length(handles.SelectedListboxContents));
    if get(handles.PlotLPF_checkbox, 'Value') || LoopGainCheckbox
        if handles.Control.Pitch.Scheduled
            for i = 1:length(handles.SelectedListboxContents)
                selIndex = findnearest(str2double(handles.SelectedListboxContents{i}), handles.Control.Pitch.ScheduledPitchAngles*180/pi);
                if handles.Control.Pitch.LowPassOrder == 1
                    Controller(1,i) = Controller(1,i)...
                        *tf(handles.Control.Pitch.LowPassCutOffFreqGS(selIndex),...
                            [1 handles.Control.Pitch.LowPassCutOffFreqGS(selIndex)]);
                else
                    Controller(1,i) = Controller(1,i)...
                        *tf(handles.Control.Pitch.LowPassCutOffFreqGS(selIndex)^2,...
                            [1 2/sqrt(2)*handles.Control.Pitch.LowPassCutOffFreqGS(selIndex) handles.Control.Pitch.LowPassCutOffFreqGS(selIndex)^2]);
                end
            end
        else
            if handles.Control.Pitch.LowPassOrder == 1
                Controller(1,:) = Controller(1,:)...
                    *tf(handles.Control.Pitch.LowPassCutOffFreq,[1 handles.Control.Pitch.LowPassCutOffFreq]);
            else
                Controller(1,:) = Controller(1,:)...
                    *tf(handles.Control.Pitch.LowPassCutOffFreq^2,...
                    [1 2/sqrt(2)*handles.Control.Pitch.LowPassCutOffFreq handles.Control.Pitch.LowPassCutOffFreq^2]);
            end
        end
    end
    if get(handles.PlotPI_checkbox, 'Value') || LoopGainCheckbox
        if handles.Control.Pitch.Scheduled
            for i = 1:length(handles.SelectedListboxContents)
                selIndex = findnearest(str2double(handles.SelectedListboxContents{i}), handles.Control.Pitch.ScheduledPitchAngles*180/pi);
                if all([handles.Control.Pitch.KpGS(selIndex) handles.Control.Pitch.KiGS(selIndex)] == 0)
                    Controller(1,i) = Controller(1,i);
                else
                    Controller(1,i) = Controller(1,i)...
                        *tf([handles.Control.Pitch.KpGS(selIndex) handles.Control.Pitch.KiGS(selIndex)], [1 0]);
                end
            end
        else
            if all([handles.Control.Pitch.Kp handles.Control.Pitch.Ki] == 0)
                Controller(1,i) = Controller(1,i);
            else
                Controller(1,:) = Controller(1,:)...
                    *tf([handles.Control.Pitch.Kp handles.Control.Pitch.Ki], [1 0]);
            end
        end
    end
    if get(handles.PlotNotch_checkbox, 'Value') || LoopGainCheckbox
        if handles.Control.Pitch.Scheduled
            for i = 1:length(handles.SelectedListboxContents)
                selIndex = findnearest(str2double(handles.SelectedListboxContents{i}), handles.Control.Pitch.ScheduledPitchAngles*180/pi);
                if any([handles.Control.Pitch.Notch_beta1GS(selIndex) handles.Control.Pitch.Notch_beta2GS(selIndex) handles.Control.Pitch.Notch_wnGS(selIndex)] == 0) ...
                   && any([handles.Control.Pitch.Notch2_beta1GS(selIndex) handles.Control.Pitch.Notch2_beta2GS(selIndex) handles.Control.Pitch.Notch2_wnGS(selIndex)] == 0)
                    Controller(1,i) = Controller(1,i);
                else
                    Controller(1,i) = Controller(1,i)...
                        *tf([1 2*handles.Control.Pitch.Notch_beta1GS(selIndex)*handles.Control.Pitch.Notch_wnGS(selIndex) handles.Control.Pitch.Notch_wnGS(selIndex)^2],...
                            [1 2*handles.Control.Pitch.Notch_beta2GS(selIndex)*handles.Control.Pitch.Notch_wnGS(selIndex) handles.Control.Pitch.Notch_wnGS(selIndex)^2])...
                        *tf([1 2*handles.Control.Pitch.Notch2_beta1GS(selIndex)*handles.Control.Pitch.Notch2_wnGS(selIndex) handles.Control.Pitch.Notch2_wnGS(selIndex)^2],...
                            [1 2*handles.Control.Pitch.Notch2_beta2GS(selIndex)*handles.Control.Pitch.Notch2_wnGS(selIndex) handles.Control.Pitch.Notch2_wnGS(selIndex)^2]);
                end
            end
        else
            if any([handles.Control.Pitch.Notch_beta1 handles.Control.Pitch.Notch_beta2 handles.Control.Pitch.Notch_wn] == 0)...
               && any([handles.Control.Pitch.Notch2_beta1 handles.Control.Pitch.Notch2_beta2 handles.Control.Pitch.Notch2_wn] == 0)
                Controller(1,:) = Controller(1,:);
            else
                Controller(1,:) = Controller(1,:)...
                    *tf([1 2*handles.Control.Pitch.Notch_beta1*handles.Control.Pitch.Notch_wn handles.Control.Pitch.Notch_wn^2], ...
                        [1 2*handles.Control.Pitch.Notch_beta2*handles.Control.Pitch.Notch_wn handles.Control.Pitch.Notch_wn^2])...
                    *tf([1 2*handles.Control.Pitch.Notch2_beta1*handles.Control.Pitch.Notch2_wn handles.Control.Pitch.Notch2_wn^2], ...
                        [1 2*handles.Control.Pitch.Notch2_beta2*handles.Control.Pitch.Notch2_wn handles.Control.Pitch.Notch2_wn^2]);
            end
        end
    end
    
function LPFCutOff_text_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.LowPassCutOffFreq))
end
handles.Control.Pitch.LowPassCutOffFreq = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LPFCutOff_text_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Table_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in SecondOrderLPF_radio.
function SecondOrderLPF_radio_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Control.Pitch.LowPassOrder = 2;
else
    handles.Control.Pitch.LowPassOrder = 1;
end
guidata(hObject, handles);


% --- Executes on button press in FirstOrderLPF_radio.
function FirstOrderLPF_radio_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Control.Pitch.LowPassOrder = 1;
else
    handles.Control.Pitch.LowPassOrder = 2;
end
guidata(hObject, handles);

function handles = UpdateHandlesWithTableData(handles)
    % Get geometry from table
    Table = get(handles.Table, 'Data');
    
    % Empty gain vector
    ScheduledPitchAngles = nan(handles.TableSize,1);
    KiGS = nan(handles.TableSize,1);
    KpGS = nan(handles.TableSize,1);
    Notch_beta1GS = nan(handles.TableSize,1);
    Notch_beta2GS = nan(handles.TableSize,1);
    Notch_wnGS = nan(handles.TableSize,1);
    % Second notch filter start
    Notch2_beta1GS = nan(handles.TableSize,1);
    Notch2_beta2GS = nan(handles.TableSize,1);
    Notch2_wnGS = nan(handles.TableSize,1);
    % Second notch filter end
    LowPassCutOffFreqGS = nan(handles.TableSize,1);

    % Find invalid cells
    for i = 1:size(Table,1)
        for j = 1:size(Table,2)
                invalid(i,j) = ...
                    isempty(cell2mat(Table(i,j))) + ...
                    sum(isnan(cell2mat(Table(i,j))));
        end
    end

    % Extract geometry from table
    for i = 1:size(Table,1)
        if invalid(i,1)
            ScheduledPitchAngles(i) = 0;
        else
            ScheduledPitchAngles(i) = Table{i,1};
        end
        if invalid(i,2)
            KpGS(i) = 0;
        else
            KpGS(i) = Table{i,2};
        end
        if invalid(i,3)
            KiGS(i) = 0;
        else
            KiGS(i) = Table{i,3};
        end
        if invalid(i,4)
            Notch_beta1GS(i) = 0;
        else
            Notch_beta1GS(i) = Table{i,4};
        end
        if invalid(i,5)
            Notch_beta2GS(i) = 0;
        else
            Notch_beta2GS(i) = Table{i,5};
        end
        if invalid(i,6)
            Notch_wnGS(i) = 0;
        else
            Notch_wnGS(i) = Table{i,6};
        end
        if invalid(i,7)
            Notch2_beta1GS(i) = 0;
        else
            Notch2_beta1GS(i) = Table{i,7};
        end
        if invalid(i,8)
            Notch2_beta2GS(i) = 0;
        else
            Notch2_beta2GS(i) = Table{i,8};
        end
        if invalid(i,9)
            Notch2_wnGS(i) = 0;
        else
            Notch2_wnGS(i) = Table{i,9};
        end
        if invalid(i,10)
            LowPassCutOffFreqGS(i) = 0;
        else
            LowPassCutOffFreqGS(i) = Table{i,10};
        end
    end
    
    handles.Control.Pitch.ScheduledPitchAngles = ScheduledPitchAngles(~isnan(ScheduledPitchAngles))*pi/180;
    handles.Control.Pitch.KpGS = KpGS(~isnan(KpGS));
    handles.Control.Pitch.KiGS = KiGS(~isnan(KiGS));
    handles.Control.Pitch.Notch_beta1GS = Notch_beta1GS(~isnan(Notch_beta1GS));
    handles.Control.Pitch.Notch_beta2GS = Notch_beta2GS(~isnan(Notch_beta2GS));
    handles.Control.Pitch.Notch_wnGS = Notch_wnGS(~isnan(Notch_wnGS));
    % Second notch filter start
    handles.Control.Pitch.Notch2_beta1GS = Notch2_beta1GS(~isnan(Notch2_beta1GS));
    handles.Control.Pitch.Notch2_beta2GS = Notch2_beta2GS(~isnan(Notch2_beta2GS));
    handles.Control.Pitch.Notch2_wnGS = Notch2_wnGS(~isnan(Notch2_wnGS));
    % Second notch filter end
    handles.Control.Pitch.LowPassCutOffFreqGS = LowPassCutOffFreqGS(~isnan(LowPassCutOffFreqGS));


function assertControllerGains(handles)
    if get(handles.GainScheduled,'Value')
        if any(handles.Control.Pitch.KpGS > 0) || any(handles.Control.Pitch.KiGS > 0)
            questdlg('Warning: Proportional and integral gains must all be negative','Warning','OK','OK');
        end
    
        if (any(handles.Control.Pitch.Notch_beta1GS < 0) || any(handles.Control.Pitch.Notch_beta2GS < 0) || any(handles.Control.Pitch.Notch_wnGS < 0)) ...
            || (any(handles.Control.Pitch.Notch2_beta1GS < 0) || any(handles.Control.Pitch.Notch2_beta2GS < 0) || any(handles.Control.Pitch.Notch2_wnGS < 0))
            questdlg('Warning: Notch parameters must all be positive','Warning','OK','OK');
        end
    
        if any(handles.Control.Pitch.LowPassCutOffFreqGS < 0)
            questdlg('Warning: Low-pass filter cut-off frequencies must all be positive','Warning','OK','OK');
        end
    else
        if handles.Control.Pitch.Kp > 0 || handles.Control.Pitch.Ki > 0
            questdlg('Warning: Proportional and integral gains must be negative','Warning','OK','OK');
        end
        
        if (any(handles.Control.Pitch.Notch_beta1 < 0) || any(handles.Control.Pitch.Notch_beta2 < 0) || any(handles.Control.Pitch.Notch_wn < 0))...
            || (any(handles.Control.Pitch.Notch2_beta1 < 0) || any(handles.Control.Pitch.Notch2_beta2 < 0) || any(handles.Control.Pitch.Notch2_wn < 0))
            questdlg('Warning: Notch parameters must be positive','Warning','OK','OK');
        end
    
        if any(handles.Control.Pitch.LowPassCutOffFreq < 0)
            questdlg('Warning: Low-pass filter cut-off frequency must be positive','Warning','OK','OK');
        end
    end
    
% --- Executes when entered data in editable cell(s) in Table.
function Table_CellEditCallback(hObject, eventdata, handles)
handles = UpdateHandlesWithTableData(handles);
guidata(hObject, handles);


% --- First Notch Filter
function Constant_Notch_B1_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Notch_beta1))
end
handles.Control.Pitch.Notch_beta1 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Notch_B1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Constant_Notch_B2_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Notch_beta2))
end
handles.Control.Pitch.Notch_beta2 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Notch_B2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Constant_Notch_wn_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Notch_wn))
end
handles.Control.Pitch.Notch_wn = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Notch_wn_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Second Notch Filter
function Constant_Notch2_B1_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Notch2_beta1))
end
handles.Control.Pitch.Notch2_beta1 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Notch2_B1_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Constant_Notch2_B2_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Notch2_beta2))
end
handles.Control.Pitch.Notch2_beta2 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Notch2_B2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Constant_Notch2_wn_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Notch2_wn))
end
handles.Control.Pitch.Notch2_wn = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Notch2_wn_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function Window_slider_Callback(hObject, eventdata, handles)
xpos = get(hObject, 'Value');
sliderpos = get(hObject, 'Position');
windowpos = get(handles.Window_group, 'Position');
windowpos(1) = xpos*(sliderpos(3)-windowpos(3));
set(handles.Window_group, 'Position', windowpos);

% --- Executes during object creation, after setting all properties.
function Window_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in ExportPlotData_pushbutton.
function ExportPlotData_pushbutton_Callback(hObject, eventdata, handles)
    BodePlot(handles, false, true)

% --- Executes on button press in RefreshPlot_pushbutton.
function RefreshPlot_pushbutton_Callback(hObject, eventdata, handles)
    BodePlot(handles, false, false)
