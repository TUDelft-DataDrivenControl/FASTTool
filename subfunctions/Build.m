%% Initialization code
function varargout = Build(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Build_OpeningFcn, ...
                   'gui_OutputFcn',  @Build_OutputFcn, ...
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
function Build_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Tower_Height = varargin{1};
handles.Tower_OuterDiameter = varargin{2};
handles.Tower_Mass = varargin{3};
handles.Tower_EI = varargin{4};
handles.Tower_BottomThickness = varargin{5};
handles.Tower_TopThickness = varargin{6};
handles.Blade_Radius = varargin{7};
handles.Blade_Mass = varargin{8};
handles.Blade_EIflap = varargin{9};
handles.Blade_EIedge = varargin{10};
handles.Blade_Twist = varargin{11};
handles.Blade_Chord = varargin{12};
handles.Blade_PitchAxis = varargin{13};
handles.Blade_NFoil = varargin{14};
handles.Blade_Cone = varargin{15};
handles.NBlades = varargin{16};
handles.Airfoils = varargin{17};
handles.AoA = varargin{18};
handles.CL = varargin{19};
handles.CD = varargin{20};
handles.Nacelle = varargin{21};
handles.Nacelle_Length = varargin{22};
handles.Nacelle_Diameter = varargin{23};
handles.Nacelle_Mass = varargin{24};
handles.Hub = varargin{25};
handles.Hub_Height = varargin{26};
handles.Hub_Overhang = varargin{27};
handles.Hub_Mass = varargin{28};
handles.Shaft_Tilt = varargin{29};
handles.Nosecone_Length = varargin{30};
handles.WindSpeed_Cutin = varargin{31};
handles.WindSpeed_Cutout = varargin{32};
handles.Rated_Power = varargin{33};
handles.Rated_TipSpeedRatio = varargin{34};
handles.Blade_PitchOffset = varargin{35};
handles.Generator_Efficiency = varargin{36};
handles.Generator_Torque = varargin{37};
handles.Generator_Speed = varargin{38};
handles.Generator_K = varargin{39};
handles.HSS_Inertia = varargin{40};
handles.Gearbox_Efficiency = varargin{41};
handles.Gearbox_Ratio = varargin{42};

% Flag for first time run
handles.FirstModalAnalysis = true;

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.Build);

%% Output function
function varargout = Build_OutputFcn(hObject, eventdata, handles) 

% Close figure
delete(hObject)

%% Close button
function Close_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume(handles.Build);

%% Modal analysis
function RunModal_Callback(hObject, eventdata, handles)

% Rotation speed
rpm = str2double(get(handles.ModesRPM_textbox, 'String'));
if isempty(rpm) || isnan(rpm)
    errordlg('Invalid rotation speed.', 'Error')
else

set(hObject, 'Enable', 'off');
    
% Get geometry from handles
Tower_Height = handles.Tower_Height;
Tower_Mass = handles.Tower_Mass;
Tower_OuterDiameter = handles.Tower_OuterDiameter;
Tower_BottomThickness = handles.Tower_BottomThickness;
Tower_TopThickness = handles.Tower_TopThickness;
Tower_WallThickness = linspace(Tower_BottomThickness,Tower_TopThickness,length(Tower_Height))';
Tower_EI = handles.Tower_EI;
Nacelle_Mass = handles.Nacelle_Mass;
Hub_Mass = handles.Hub_Mass;
Blade_Radius = handles.Blade_Radius;
Blade_Chord = handles.Blade_Chord;
Blade_Mass = handles.Blade_Mass;
Blade_EIedge = handles.Blade_EIedge;
Blade_EIflap = handles.Blade_EIflap;
Blade_Twist = handles.Blade_Twist;
Blade_NFoil = handles.Blade_NFoil;
NBlades = handles.NBlades;
Airfoils = handles.Airfoils;
    
% Only calculate tower modes for 0 rpm (first time run)
if handles.FirstModalAnalysis

% Tower structural properties
Shear_Modulus = 80.8e9;
Youngs_Modulus = 210e9;
Tower_GJ = Shear_Modulus * pi/32*(Tower_OuterDiameter.^4 - (Tower_OuterDiameter-2*Tower_WallThickness).^4);
Tower_EA = Youngs_Modulus * pi*Tower_OuterDiameter.*Tower_WallThickness;
Tower_Iner = 0.5 * Tower_Mass .* (Tower_OuterDiameter/2 - Tower_WallThickness).^2;

% Tower input file
n_secs = length(Tower_Height);
fid = fopen('subfunctions\tower_sec_props.dat', 'wt');
fprintf(fid, 'Tower section properties\n');
fprintf(fid, '%0.0f\tn_secs:\tnumber of blade or tower sections at which properties are specified (-)\n\n', n_secs);
fprintf(fid, 'sec_loc\tstr_tw\ttw_iner\tmass_den\tflp_iner\tedge_iner\tflp_stff\tedge_stff\ttor_stff\taxial_stff\tcg_offst\tsc_offst\ttc_offst\n');
fprintf(fid, '(-)\t(deg)\t(deg)\t(kg/m)\t\t(kg-m)\t\t(kg-m)\t\t(Nm^2)\t\t(Nm^2)\t\t(Nm^2)\t\t(N)\t\t(m)\t\t(m)\t\t(m)\n');
for i = 1:n_secs
    fprintf(fid, '%0.5f\t%0.0f\t%0.0f\t%7.2f\t\t%8.2f\t%8.2f\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\n', ...
        Tower_Height(i)/Tower_Height(end), ...
        0, ...
        0, ...
        Tower_Mass(i), ...
        Tower_Iner(i), ...
        Tower_Iner(i), ...
        Tower_EI(i), ...
        Tower_EI(i), ...
        Tower_GJ(i), ...
        Tower_EA(i), ...
        0, ...
        0, ...
        0);
end
fclose(fid);

% BModes tower input file
fid = fopen('subfunctions\tower.bmi', 'wt');
fprintf(fid, '======================   BModes v1.03 Main Input File  ==================\n');
fprintf(fid, 'Sample non-uniform blade (output is space-delimited)\n\n');
fprintf(fid, '--------- General parameters ---------------------------------------------------------------------\n');
fprintf(fid, 'False    Echo        Echo input file contents to *.echo file if true.\n');
fprintf(fid, '2       beam_type   1: blade, 2: tower (-)\n');
fprintf(fid, '%4.1f    romg:       rotor speed, automatically set to zero for tower modal analysis (rpm)\n', rpm);
fprintf(fid, '1.0      romg_mult:  rotor speed muliplicative factor (-)\n');
fprintf(fid, '%4.2f    radius:     rotor tip radius measured along coned blade axis OR tower height (m)\n', Tower_Height(end));
fprintf(fid, '%4.2f    hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\n', Tower_Height(1));
fprintf(fid, '%4.2f    precone:    built-in precone angle, automatically set to zero for a tower (deg)\n', handles.Blade_Cone);
fprintf(fid, '%4.3f    bl_thp:     blade pitch setting, automatically set to zero for a tower (deg)\n', handles.Blade_PitchOffset);
fprintf(fid, '1         hub_conn:   hub-to-blade connection [1: cantilevered; other options not yet available] (-)\n');
fprintf(fid, '20        modepr:     number of modes to be printed (-)\n');
fprintf(fid, 't         TabDelim    (true: tab-delimited output tables; false: space-delimited tables)\n');
fprintf(fid, 't         mid_node_tw  (true: output twist at mid-node of elements; false: no mid-node outputs)\n\n');
fprintf(fid, '--------- Blade-tip or tower-top mass properties --------------------------------------------\n');
fprintf(fid, '%4.2f    tip_mass    blade-tip or tower-top mass (kg)\n', Nacelle_Mass + Hub_Mass + NBlades*trapz(Blade_Radius,Blade_Mass));
fprintf(fid, '0.        cm_loc      tip-mass c.m. offset from the blade axis measured along the tip section y reference axis (m)\n');
fprintf(fid, '0.        ixx_tip     blade lag mass moment of inertia about the tip-section x reference axis (kg-m^2)\n');
fprintf(fid, '0.        iyy_tip     blade flap mass moment of inertia about the tip-section y reference axis (kg-m^2)\n');
fprintf(fid, '0.        izz_tip     torsion mass moment of inertia about the tip-section z reference axis (kg-m^2)\n');
fprintf(fid, '0.        ixy_tip     cross product of inertia about x and y reference axes(kg-m^2)\n');
fprintf(fid, '0.        izx_tip     cross product of inertia about z and x reference axes(kg-m^2)\n');
fprintf(fid, '0.        iyz_tip     cross product of inertia about y and z reference axes(kg-m^2)\n\n');
fprintf(fid, '--------- Distributed-property identifiers --------------------------------------------------------\n');
fprintf(fid, '1         id_mat:     material_type [1: isotropic; non-isotropic composites option not yet available]\n');
fprintf(fid, '%s sec_props_file   name of beam section properties file (-)\n\n', 'subfunctions\tower_sec_props.dat');
fprintf(fid, 'Property scaling factors..............................\n');
fprintf(fid, '1.0       sec_mass_mult:   mass density multiplier (-)\n');
fprintf(fid, '1.0       flp_iner_mult:   blade flap or tower f-a inertia multiplier (-)\n');
fprintf(fid, '1.0       lag_iner_mult:   blade lag or tower s-s inertia multiplier (-)\n');
fprintf(fid, '1.0       flp_stff_mult:   blade flap or tower f-a bending stiffness multiplier (-)\n');
fprintf(fid, '1.0       edge_stff_mult:  blade lag or tower s-s bending stiffness multiplier (-)\n');
fprintf(fid, '1.0       tor_stff_mult:   torsion stiffness multiplier (-)\n');
fprintf(fid, '1.0       axial_stff_mult: axial stiffness multiplier (-)\n');
fprintf(fid, '1.0       cg_offst_mult:   cg offset multiplier (-)\n');
fprintf(fid, '1.0       sc_offst_mult:   shear center multiplier (-)\n');
fprintf(fid, '1.0       tc_offst_mult:   tension center multiplier (-)\n\n');
fprintf(fid, '--------- Finite element discretization --------------------------------------------------\n');
fprintf(fid, '20        nselt:     no of blade or tower elements (-)\n');
fprintf(fid, 'Distance of element boundary nodes from blade or flexible-tower root (normalized wrt blade or tower length), el_loc()\n');
fprintf(fid, '0.0	0.05	0.1	0.15	0.2	0.25	0.3	0.35	0.4	0.45	0.5	0.55	0.6	0.65	0.7	0.75	0.8	0.85	0.9	0.95	1.0\n\n');
fprintf(fid, '--------- Properties of tension wires suporting the tower --------------------------------\n');
fprintf(fid, '0         n_attachments: no of wire-attachment locations on tower, maxm allowable is 2; 0: no tension-wire support (-)\n');
fprintf(fid, '3 3       n_wires:       no of wires attached at each location (must be 3 or higher) (-)\n');
fprintf(fid, '6 9       node_attach:   node numbers of attacments location (node number must be more than 1 and less than nselt+2) (-)\n');
fprintf(fid, '0.e0 0.e0 wire_stfness:  wire spring constant in each set (N/m)\n');
fprintf(fid, '0. 0.     th_wire:       angle of tension wires wrt the tower axis at each attachment point (deg)\n');
fclose(fid);

% Run BModes
system('subfunctions\bmodes subfunctions\tower.bmi');
    
% Read tower output
fileID = fopen('subfunctions\tower.out','r');
textscan(fileID, '%[^\n\r]', 5, 'ReturnOnError', false);
data = textscan(fileID, '%s%s%s%s%s%s%[^\n\r]', 'Delimiter', '\t', 'ReturnOnError', false);
fclose(fileID);
n = (length(data{1}) - 1)/20;
mode = zeros(20,1);
freq = zeros(20,1);
y_SS = zeros(20,n-2);
y_FA = zeros(20,n-2);
y_tw = zeros(20,n-2);
for i = 1:20
    header = cell2mat(data{1}((i-1)*n+1));
    freq(i) = str2double(header(end-15:end-4));
    for j = 1:n-2
        y_SS(i,j) = str2double(data{2}((i-1)*n+2+j));
        y_FA(i,j) = str2double(data{4}((i-1)*n+2+j));
        y_tw(i,j) = pi/180*str2double(data{6}((i-1)*n+2+j));
    end
    if max(abs(y_SS(i,:))) > max(abs(y_FA(i,:))) && max(abs(y_SS(i,:))) > max(abs(y_tw(i,:)))
        mode(i) = 1;
    elseif max(abs(y_FA(i,:))) > max(abs(y_SS(i,:))) && max(abs(y_FA(i,:))) > max(abs(y_tw(i,:)))
        mode(i) = 2;
    else
        mode(i) = 3;
    end
end
y_SS = y_SS(:,1:2:end);
y_FA = y_FA(:,1:2:end);
y_tw = y_tw(:,1:2:end);

% Tower side-side modes
i = find(mode == 1);
x = linspace(0,1,(n-3)/2+1);
Tower_SideSide1_shape = [y_FA(i(1),:); y_SS(i(1),:); y_tw(i(1),:)];
Tower_SideSide1_coeff = polyfit(x,y_SS(i(1),:),6);
Tower_SideSide1_coeff = Tower_SideSide1_coeff / sum(Tower_SideSide1_coeff(1:5));
Tower_SideSide1_freq = freq(i(1));
Tower_SideSide2_shape = [y_FA(i(2),:); y_SS(i(2),:); y_tw(i(2),:)];
Tower_SideSide2_coeff = polyfit(x,y_SS(i(2),:),6);
Tower_SideSide2_coeff = Tower_SideSide2_coeff / sum(Tower_SideSide2_coeff(1:5));
Tower_SideSide2_freq = freq(i(2));

% Tower fore-aft modes
i = find(mode == 2);
x = linspace(0,1,(n-3)/2+1);
Tower_ForeAft1_shape = [y_FA(i(1),:); y_SS(i(1),:); y_tw(i(1),:)];
Tower_ForeAft1_coeff = polyfit(x,y_FA(i(1),:),6);
Tower_ForeAft1_coeff = Tower_ForeAft1_coeff / sum(Tower_ForeAft1_coeff(1:5));
Tower_ForeAft1_freq = freq(i(1));
Tower_ForeAft2_shape = [y_FA(i(2),:); y_SS(i(2),:); y_tw(i(2),:)];
Tower_ForeAft2_coeff = polyfit(x,y_FA(i(2),:),6);
Tower_ForeAft2_coeff = Tower_ForeAft2_coeff / sum(Tower_ForeAft2_coeff(1:5));
Tower_ForeAft2_freq = freq(i(2));

end

% Blade Thickness
Blade_Thickness = zeros(size(Blade_Radius));
for i = 1:length(Blade_NFoil)
    t_u = max(Airfoils{Blade_NFoil(i)}(2,1:200) * Blade_Chord(i));
    t_l = min(Airfoils{Blade_NFoil(i)}(2,200:end) * Blade_Chord(i));
    Blade_Thickness(i) = t_u - t_l;
end
s = 2.*round((length(Blade_NFoil)/5+1)/2)-1;
Blade_Thickness = conv(Blade_Thickness(:),s(:)/sum(s),'same')./Blade_Thickness(:);

% Offsets
i = find(Blade_NFoil > 2);
x = [Blade_Radius(1), Blade_Radius(i(1)), Blade_Radius(end)];
cg = [0, 0.2, 0.2];
sc = [0, -0.03, 0.1];

% Factors for inertia and torsional stiffness per airfoil
Cflap = [0.446; 0.260; 0.147; 0.035; 0.027; 0.014; 0.010; 0.004];
Cedge = [0.034; 0.028; 0.022; 0.019; 0.017; 0.014; 0.015; 0.018];
Ctor =  [0.170; 0.136; 0.094; 0.020; 0.021; 0.025; 0.025; 0.022];
Ccgo =  [0.010; 0.018; 0.030; 0.060; 0.047; 0.037; 0.045; 0.060];

% Blade structural properties
Blade_FlapIner = Cflap(Blade_NFoil).*(Blade_Mass.*Blade_Thickness.^2);
Blade_EdgeIner = Cedge(Blade_NFoil).*(Blade_Mass.*Blade_Chord.^2);
Blade_GJ = Ctor(Blade_NFoil).*(Blade_EIflap+Blade_EIedge);
Blade_EA = 1.3e7 * Blade_Mass;
Blade_cg = interp1(x,cg,Blade_Radius) .* Blade_Chord;
Blade_sc = interp1(x,sc,Blade_Radius) .* Blade_Chord;

% Blade aerodynamic properties
Blade_ac = NaN(size(Blade_Radius));
Blade_ac(Blade_NFoil == 1) = 0.25;
Blade_ac(Blade_NFoil >= 4) = 0.125;
Blade_ac(isnan(Blade_ac)) = interp1(...
    Blade_Radius(~isnan(Blade_ac)), ...
    Blade_ac(~isnan(Blade_ac)), ...
    Blade_Radius(isnan(Blade_ac)));
Blade_eo = Ccgo(Blade_NFoil) .* Blade_Chord;

% Blade input file
n_secs = length(Blade_Radius);
fid = fopen('subfunctions\blade_sec_props.dat', 'wt');
fprintf(fid, 'Blade section properties\n');
fprintf(fid, '%0.0f\tn_secs:\tnumber of blade or tower sections at which properties are specified (-)\n\n', n_secs);
fprintf(fid, 'sec_loc\tstr_tw\ttw_iner\tmass_den\tflp_iner\tedge_iner\tflp_stff\tedge_stff\ttor_stff\taxial_stff\tcg_offst\tsc_offst\ttc_offst\n');
fprintf(fid, '(-)\t(deg)\t(deg)\t(kg/m)\t\t(kg-m)\t\t(kg-m)\t\t(Nm^2)\t\t(Nm^2)\t\t(Nm^2)\t\t(N)\t\t(m)\t\t(m)\t\t(m)\n');
for i = 1:n_secs
    fprintf(fid, '%0.5f\t%7.3f\t%7.3f\t%7.2f\t\t%8.2f\t%8.2f\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\n', ...
        (Blade_Radius(i)-Blade_Radius(1))/(Blade_Radius(end)-Blade_Radius(1)), ...
        Blade_Twist(i), ...
        Blade_Twist(i), ...
        Blade_Mass(i), ...
        Blade_FlapIner(i), ...
        Blade_EdgeIner(i), ...
        Blade_EIflap(i), ...
        Blade_EIedge(i), ...
        Blade_GJ(i), ...
        Blade_EA(i), ...
        Blade_cg(i), ...
        Blade_sc(i), ...
        0);
end
fclose(fid);

% BModes blade input file  
fid = fopen('subfunctions\blade.bmi', 'wt');
fprintf(fid, '======================   BModes v1.03 Main Input File  ==================\n');
fprintf(fid, 'Sample non-uniform blade (output is space-delimited)\n\n');
fprintf(fid, '--------- General parameters ---------------------------------------------------------------------\n');
fprintf(fid, 'False    Echo        Echo input file contents to *.echo file if true.\n');
fprintf(fid, '1       beam_type   1: blade, 2: tower (-)\n', i);
fprintf(fid, '%4.1f    romg:       rotor speed, automatically set to zero for tower modal analysis (rpm)\n', rpm);
fprintf(fid, '1.0      romg_mult:  rotor speed muliplicative factor (-)\n');
fprintf(fid, '%4.2f    radius:     rotor tip radius measured along coned blade axis OR tower height (m)\n', Blade_Radius(end));
fprintf(fid, '%4.2f    hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\n', Blade_Radius(1));
fprintf(fid, '%4.2f    precone:    built-in precone angle, automatically set to zero for a tower (deg)\n', handles.Blade_Cone);
fprintf(fid, '%4.2f    bl_thp:     blade pitch setting, automatically set to zero for a tower (deg)\n', handles.Blade_PitchOffset);
fprintf(fid, '1         hub_conn:   hub-to-blade connection [1: cantilevered; other options not yet available] (-)\n');
fprintf(fid, '20        modepr:     number of modes to be printed (-)\n');
fprintf(fid, 't         TabDelim    (true: tab-delimited output tables; false: space-delimited tables)\n');
fprintf(fid, 't         mid_node_tw  (true: output twist at mid-node of elements; false: no mid-node outputs)\n\n');
fprintf(fid, '--------- Blade-tip or tower-top mass properties --------------------------------------------\n');
fprintf(fid, '%4.2f    tip_mass    blade-tip or tower-top mass (kg)\n', 0);
fprintf(fid, '0.        cm_loc      tip-mass c.m. offset from the blade axis measured along the tip section y reference axis (m)\n');
fprintf(fid, '0.        ixx_tip     blade lag mass moment of inertia about the tip-section x reference axis (kg-m^2)\n');
fprintf(fid, '0.        iyy_tip     blade flap mass moment of inertia about the tip-section y reference axis (kg-m^2)\n');
fprintf(fid, '0.        izz_tip     torsion mass moment of inertia about the tip-section z reference axis (kg-m^2)\n');
fprintf(fid, '0.        ixy_tip     cross product of inertia about x and y reference axes(kg-m^2)\n');
fprintf(fid, '0.        izx_tip     cross product of inertia about z and x reference axes(kg-m^2)\n');
fprintf(fid, '0.        iyz_tip     cross product of inertia about y and z reference axes(kg-m^2)\n\n');
fprintf(fid, '--------- Distributed-property identifiers --------------------------------------------------------\n');
fprintf(fid, '1         id_mat:     material_type [1: isotropic; non-isotropic composites option not yet available]\n');
fprintf(fid, '%s sec_props_file   name of beam section properties file (-)\n\n', 'subfunctions\blade_sec_props.dat');
fprintf(fid, 'Property scaling factors..............................\n');
fprintf(fid, '1.0       sec_mass_mult:   mass density multiplier (-)\n');
fprintf(fid, '1.0       flp_iner_mult:   blade flap or tower f-a inertia multiplier (-)\n');
fprintf(fid, '1.0       lag_iner_mult:   blade lag or tower s-s inertia multiplier (-)\n');
fprintf(fid, '1.0       flp_stff_mult:   blade flap or tower f-a bending stiffness multiplier (-)\n');
fprintf(fid, '1.0       edge_stff_mult:  blade lag or tower s-s bending stiffness multiplier (-)\n');
fprintf(fid, '1.0       tor_stff_mult:   torsion stiffness multiplier (-)\n');
fprintf(fid, '1.0       axial_stff_mult: axial stiffness multiplier (-)\n');
fprintf(fid, '1.0       cg_offst_mult:   cg offset multiplier (-)\n');
fprintf(fid, '1.0       sc_offst_mult:   shear center multiplier (-)\n');
fprintf(fid, '1.0       tc_offst_mult:   tension center multiplier (-)\n\n');
fprintf(fid, '--------- Finite element discretization --------------------------------------------------\n');
fprintf(fid, '20        nselt:     no of blade or tower elements (-)\n');
fprintf(fid, 'Distance of element boundary nodes from blade or flexible-tower root (normalized wrt blade or tower length), el_loc()\n');
fprintf(fid, '0.0	0.05	0.1	0.15	0.2	0.25	0.3	0.35	0.4	0.45	0.5	0.55	0.6	0.65	0.7	0.75	0.8	0.85	0.9	0.95	1.0\n\n');
fprintf(fid, '--------- Properties of tension wires suporting the tower --------------------------------\n');
fprintf(fid, '0         n_attachments: no of wire-attachment locations on tower, maxm allowable is 2; 0: no tension-wire support (-)\n');
fprintf(fid, '3 3       n_wires:       no of wires attached at each location (must be 3 or higher) (-)\n');
fprintf(fid, '6 9       node_attach:   node numbers of attacments location (node number must be more than 1 and less than nselt+2) (-)\n');
fprintf(fid, '0.e0 0.e0 wire_stfness:  wire spring constant in each set (N/m)\n');
fprintf(fid, '0. 0.     th_wire:       angle of tension wires wrt the tower axis at each attachment point (deg)\n');
fclose(fid);

% Run BModes
system('subfunctions\bmodes subfunctions\blade.bmi');

% Read blade output
fileID = fopen('subfunctions\blade.out','r');
textscan(fileID, '%[^\n\r]', 5, 'ReturnOnError', false);
data = textscan(fileID, '%s%s%s%s%s%s%[^\n\r]', 'Delimiter', '\t', 'ReturnOnError', false);
fclose(fileID);
n = (length(data{1}) - 1)/20;
mode = zeros(20,1);
freq = zeros(20,1);
y_flap = zeros(20,n-2);
y_edge = zeros(20,n-2);
y_tw = zeros(20,n-2);
for i = 1:20
    header = cell2mat(data{1}((i-1)*n+1));
    freq(i) = str2double(header(end-15:end-4));
    for j = 1:n-2
        y_flap(i,j) = str2double(data{2}((i-1)*n+2+j));
        y_edge(i,j) = str2double(data{4}((i-1)*n+2+j));
        y_tw(i,j) = pi/180*str2double(data{6}((i-1)*n+2+j));
    end
    if max(abs(y_flap(i,:))) > max(abs(y_edge(i,:))) && max(abs(y_flap(i,:))) > max(abs(y_tw(i,:)))
        mode(i) = 1;
    elseif max(abs(y_edge(i,:))) > max(abs(y_flap(i,:))) && max(abs(y_edge(i,:))) > max(abs(y_tw(i,:)))
        mode(i) = 2;
    else
        mode(i) = 3;
    end
end
y_flap = y_flap(:,1:2:end);
y_edge = y_edge(:,1:2:end);
y_tw = y_tw(:,1:2:end);

% Blade flapwise modes
i = find(mode == 1);
x = linspace(0,1,(n-3)/2+1);
Blade_Flap1_shape = [y_flap(i(1),:); y_edge(i(1),:); y_tw(i(1),:)];
Blade_Flap1_coeff = polyfit(x,y_flap(i(1),:),6);
Blade_Flap1_coeff = Blade_Flap1_coeff / sum(Blade_Flap1_coeff(1:5));
Blade_Flap1_freq = freq(i(1));
Blade_Flap2_shape = [y_flap(i(2),:); y_edge(i(2),:); y_tw(i(2),:)];
Blade_Flap2_coeff = polyfit(x,y_flap(i(2),:),6);
Blade_Flap2_coeff = Blade_Flap2_coeff / sum(Blade_Flap2_coeff(1:5));
Blade_Flap2_freq = freq(i(2));

% Blade edgewise modes
i = find(mode == 2);
x = linspace(0,1,(n-3)/2+1);
Blade_Edge1_shape = [y_flap(i(1),:); y_edge(i(1),:); y_tw(i(1),:)];
Blade_Edge1_coeff = polyfit(x,y_edge(i(1),:),6);
Blade_Edge1_coeff = Blade_Edge1_coeff / sum(Blade_Edge1_coeff(1:5));
Blade_Edge1_freq = freq(i(1));
Blade_Edge2_shape = [y_flap(i(2),:); y_edge(i(2),:); y_tw(i(2),:)];
Blade_Edge2_coeff = polyfit(x,y_edge(i(2),:),6);
Blade_Edge2_coeff = Blade_Edge2_coeff / sum(Blade_Edge2_coeff(1:5));
Blade_Edge2_freq = freq(i(2));

% Update fields
set(handles.Blade_Flap1_textbox, 'String', num2str(Blade_Flap1_freq, '%0.2f'));
set(handles.Blade_Flap2_textbox, 'String', num2str(Blade_Flap2_freq, '%0.2f'));
set(handles.Blade_Edge1_textbox, 'String', num2str(Blade_Edge1_freq, '%0.2f'));
set(handles.Blade_Edge2_textbox, 'String', num2str(Blade_Edge2_freq, '%0.2f'));

% Update blade modes
handles.Blade_FlapIner = Blade_FlapIner;
handles.Blade_EdgeIner = Blade_EdgeIner;
handles.Blade_GJ = Blade_GJ;
handles.Blade_EA = Blade_EA;
handles.Blade_cg = Blade_eo;
handles.Blade_ac = Blade_ac;
handles.Blade_Flap1_shape = Blade_Flap1_shape;
handles.Blade_Flap2_shape = Blade_Flap2_shape;
handles.Blade_Flap1_freq = Blade_Flap1_freq;
handles.Blade_Flap2_freq = Blade_Flap2_freq;
handles.Blade_Edge1_shape = Blade_Edge1_shape;
handles.Blade_Edge2_shape = Blade_Edge2_shape;
handles.Blade_Edge1_freq = Blade_Edge1_freq;
handles.Blade_Edge2_freq = Blade_Edge2_freq;

% If this is a first modal analysis
if handles.FirstModalAnalysis
    
    % Update fields
    set(handles.Tower_ForeAft1_textbox, 'String', num2str(Tower_ForeAft1_freq, '%0.2f'));
    set(handles.Tower_ForeAft2_textbox, 'String', num2str(Tower_ForeAft2_freq, '%0.2f'));
    set(handles.Tower_SideSide1_textbox, 'String', num2str(Tower_SideSide1_freq, '%0.2f'));
    set(handles.Tower_SideSide2_textbox, 'String', num2str(Tower_SideSide2_freq, '%0.2f'));
    
    % Update tower modes
    handles.Tower_GJ = Tower_GJ;
    handles.Tower_EA = Tower_EA;
    handles.Tower_Iner = Tower_Iner;
    handles.Tower_ForeAft1_shape = Tower_ForeAft1_shape;
    handles.Tower_ForeAft2_shape = Tower_ForeAft2_shape;
    handles.Tower_ForeAft1_freq = Tower_ForeAft1_freq;
    handles.Tower_ForeAft2_freq = Tower_ForeAft2_freq;
    handles.Tower_SideSide1_shape = Tower_SideSide1_shape;
    handles.Tower_SideSide2_shape = Tower_SideSide2_shape;
    handles.Tower_SideSide1_freq = Tower_SideSide1_freq;
    handles.Tower_SideSide2_freq = Tower_SideSide2_freq;
    
    % Store mode polynomials for 0 rpm
    handles.Tower_ForeAft1_coeff = Tower_ForeAft1_coeff;
    handles.Tower_ForeAft2_coeff = Tower_ForeAft2_coeff;
    handles.Tower_SideSide1_coeff = Tower_SideSide1_coeff;
    handles.Tower_SideSide2_coeff = Tower_SideSide2_coeff;
    handles.Blade_Flap1_coeff = Blade_Flap1_coeff;
    handles.Blade_Flap2_coeff = Blade_Flap2_coeff;
    handles.Blade_Edge1_coeff = Blade_Edge1_coeff;
    handles.Blade_Edge2_coeff = Blade_Edge2_coeff;
    
    % Enable buttons
    set(handles.SetFASTname, 'Enable', 'on');
    set(handles.ModesRPM_textbox, 'Enable', 'on');
    set(handles.Plot_Tower_ForeAft1, 'Enable', 'on');
    set(handles.Plot_Tower_ForeAft2, 'Enable', 'on');
    set(handles.Plot_Tower_SideSide1, 'Enable', 'on');
    set(handles.Plot_Tower_SideSide2, 'Enable', 'on');
    set(handles.Plot_Blade_Flap1, 'Enable', 'on');
    set(handles.Plot_Blade_Flap2, 'Enable', 'on');
    set(handles.Plot_Blade_Edge1, 'Enable', 'on');
    set(handles.Plot_Blade_Edge2, 'Enable', 'on');
    set(handles.TurbineFile_textbox, 'Enable', 'inactive');
    set(handles.TowerFile_textbox, 'Enable', 'inactive');
    set(handles.BladeFile_textbox, 'Enable', 'inactive');
    set(handles.AeroDynFile_textbox, 'Enable', 'inactive');
    set(handles.ADAMSFile_textbox, 'Enable', 'inactive');
    set(handles.LinearFile_textbox, 'Enable', 'inactive');
    
    % Flag
    handles.FirstModalAnalysis = false;
    
end

set(hObject, 'Enable', 'on');

% Update handles structure
guidata(hObject, handles);

end

%% Plot tower mode shapes
function PlotTowerMode(shape,mode,handles)

% Get geometry from handles
Tower_Height = handles.Tower_Height;
Tower_OuterDiameter = handles.Tower_OuterDiameter;

% Create fixed number of tower elements (N = 21)
N = 21;
x = linspace(Tower_Height(1), Tower_Height(end), N);
Tower_OuterDiameter = interp1(Tower_Height, Tower_OuterDiameter, x);
Tower_Height = x;

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
r = repmat(Tower_OuterDiameter(:)/2,[1,N]);
azi = repmat(linspace(0,2*pi,N),[length(Tower_OuterDiameter),1]);
x = r.*cos(azi);
y = r.*sin(azi);
z = repmat(Tower_Height(:),[1,N]);
surf(x,y,z, ...
    'FaceColor', [1, 1, 1], ...
    'EdgeColor', [0.75, 0.75, 0.75], ...
    'AmbientStrength', 0.8, ...
    'DiffuseStrength', 0.5, ...
    'SpecularStrength', 0.9, ...
    'BackFaceLighting', 'reverselit', ...
    'FaceAlpha', 0.25)

% Wall
vert = [[-1 1 1 -1 -1 1 1 -1] * Tower_OuterDiameter(1); ...
        [-1 -1 1 1 -1 -1 1 1] * Tower_OuterDiameter(1); ...
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
r = repmat(Tower_OuterDiameter(:)/2,[1,N]);
azi = repmat(linspace(0,2*pi,N),[length(Tower_OuterDiameter),1]);
x = r.*cos(azi) + repmat(Tower_Height(end)*shape(1,:)',[1,N]);
y = r.*sin(azi) + repmat(Tower_Height(end)*shape(2,:)',[1,N]);
z = repmat(Tower_Height(:),[1,N]);
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
Blade_Radius = handles.Blade_Radius;
Blade_Chord = handles.Blade_Chord;
Blade_Twist = handles.Blade_Twist;
Blade_NFoil = handles.Blade_NFoil;
Airfoils = handles.Airfoils;

% Create fixed number of blade elements (N = 21)
N = 21;
x = linspace(Blade_Radius(1), Blade_Radius(end), N);
Blade_Twist = interp1(Blade_Radius, Blade_Twist, x);
Blade_Chord = interp1(Blade_Radius, Blade_Chord, x);
Blade_NFoil = interp1(Blade_Radius, Blade_NFoil, x, 'nearest');
Blade_Radius = x;

% Smooth pitch axis
Blade_NFoil(Blade_NFoil < 1) = 1;
Blade_PitchAxis = NaN(size(Blade_Radius));
Blade_PitchAxis(Blade_NFoil == 1) = 0.5;
Blade_PitchAxis(Blade_NFoil >= 4) = 0.375;
Blade_PitchAxis(isnan(Blade_PitchAxis)) = interp1(...
    Blade_Radius(~isnan(Blade_PitchAxis)), ...
    Blade_PitchAxis(~isnan(Blade_PitchAxis)), ...
    Blade_Radius(isnan(Blade_PitchAxis)));

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
for i = 1:length(Blade_NFoil)
    x = [x; ([Airfoils{Blade_NFoil(i)}(1,1:N:end-1), Airfoils{Blade_NFoil(i)}(1,1)] - Blade_PitchAxis(i)) * Blade_Chord(i)];
end
for i = 1:length(Blade_NFoil)
    t_u = max(Airfoils{Blade_NFoil(i)}(2,1:200) * Blade_Chord(i));
    t_l = min(Airfoils{Blade_NFoil(i)}(2,200:end) * Blade_Chord(i));
    t(i) = t_u - t_l;
end
s = 2.*round((length(Blade_NFoil)/5+1)/2)-1;
t = conv(t(:),s(:)/sum(s),'same')./t(:);
for i = 1:length(Blade_NFoil)
    y = [y; -1*[Airfoils{Blade_NFoil(i)}(2,1:N:end-1), Airfoils{Blade_NFoil(i)}(2,1)] * Blade_Chord(i)*t(i)];
end
z = repmat(Blade_Radius(:),[1,size(x,2)]);
for i = 1:length(Blade_NFoil)
    t = pi/2 - Blade_Twist(i) * pi/180;
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
vert = [[-1 -1 -1 -1 0 0 0 0]/10 + Blade_Radius(1); ...
        [-1 -1 1 1 -1 -1 1 1] * Blade_Chord(1); ...
        [-1 1 1 -1 -1 1 1 -1] * Blade_Chord(1)];
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
for i = 1:length(Blade_NFoil)
    x = [x; ([Airfoils{Blade_NFoil(i)}(1,1:N:end-1), Airfoils{Blade_NFoil(i)}(1,1)] - Blade_PitchAxis(i)) * Blade_Chord(i)];
end
for i = 1:length(Blade_NFoil)
    t_u = max(Airfoils{Blade_NFoil(i)}(2,1:200) * Blade_Chord(i));
    t_l = min(Airfoils{Blade_NFoil(i)}(2,200:end) * Blade_Chord(i));
    t(i) = t_u - t_l;
end
s = 2.*round((length(Blade_NFoil)/5+1)/2)-1;
t = conv(t(:),s(:)/sum(s),'same')./t(:);
for i = 1:length(Blade_NFoil)
    y = [y; -1*[Airfoils{Blade_NFoil(i)}(2,1:N:end-1), Airfoils{Blade_NFoil(i)}(2,1)] * Blade_Chord(i)*t(i)];
end
z = repmat(Blade_Radius(:),[1,size(x,2)]);
for i = 1:length(Blade_NFoil)
    t = pi/2 - Blade_Twist(i) * 2*pi/180 - pi/180*shape(3,i);
    Rz = [cos(t),-sin(t), 0; ...
          sin(t), cos(t), 0; ...
          0,      0,      1];
    A = Rz * [x(i,:); y(i,:); z(i,:)];
    x(i,:) = A(1,:) + Blade_Radius(end)*shape(1,i);
    y(i,:) = A(2,:) + Blade_Radius(end)*shape(2,i);
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

%% Modal analysis rotation speed - textbox
function ModesRPM_textbox_Callback(hObject, eventdata, handles)
function ModesRPM_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Set file names
function SetFASTname_Callback(hObject, eventdata, handles)

% Set file name
[FileName,PathName] = uiputfile('*.fst', 'Set names for FAST input files');

if FileName

    % Update text boxes
    set(handles.TurbineFile_textbox, 'String', [FileName]);
    set(handles.TowerFile_textbox, 'String', [FileName(1:end-4), '_tower.dat']);
    set(handles.BladeFile_textbox, 'String', [FileName(1:end-4), '_blade.dat']);
    set(handles.AeroDynFile_textbox, 'String', [FileName(1:end-4), '_AeroDyn.ipt']);
    set(handles.ADAMSFile_textbox, 'String', [FileName(1:end-4), '_ADAMS.dat']);
    set(handles.LinearFile_textbox, 'String', [FileName(1:end-4), '_linearized.mat']);

    % Enable buttons
    set(handles.WriteFAST, 'Enable', 'on');
    set(handles.TurbineFile_textbox, 'Enable', 'inactive');
    set(handles.TowerFile_textbox, 'Enable', 'inactive');
    set(handles.BladeFile_textbox, 'Enable', 'inactive');
    set(handles.AeroDynFile_textbox, 'Enable', 'inactive');
    set(handles.ADAMSFile_textbox, 'Enable', 'inactive');
    set(handles.LinearFile_textbox, 'Enable', 'inactive');
    
    % Store file name
    handles.FileName = FileName(1:end-4);
    handles.PathName = PathName;
    guidata(hObject, handles);

end
function TurbineFile_textbox_Callback(hObject, eventdata, handles)
function TurbineFile_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TowerFile_textbox_Callback(hObject, eventdata, handles)
function TowerFile_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function BladeFile_textbox_Callback(hObject, eventdata, handles)
function BladeFile_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function AeroDynFile_textbox_Callback(hObject, eventdata, handles)
function AeroDynFile_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ADAMSFile_textbox_Callback(hObject, eventdata, handles)
function ADAMSFile_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LinearFile_textbox_Callback(hObject, eventdata, handles)
function LinearFile_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Write FAST input files
function WriteFAST_Callback(hObject, eventdata, handles)

% Get file name
PathName = handles.PathName;
FileName = handles.FileName;

% Prevention of some common errors:
% (1) Ensure that tower top < hub height (to prevent negative Twr2Sft)
if handles.Tower_Height(end) >= handles.Hub_Height
    handles.Tower_Height(end) = handles.Hub_Height - 0.1;
end
% (2) Blade nodes in 3-decimal precision (to prevent DRNodes not matching with RNodes in file)
r = 1e-3 * round(handles.Blade_Radius*1e3);
Blade_r = (r(1:end-1) + r(2:end))/2;
Blade_dr = r(2:end) - r(1:end-1);
Blade_NFoil = interp1(handles.Blade_Radius,handles.Blade_NFoil,Blade_r,'nearest');
Blade_Twist = interp1(handles.Blade_Radius,handles.Blade_Twist,Blade_r,'nearest');
Blade_Chord = interp1(handles.Blade_Radius,handles.Blade_Chord,Blade_r,'nearest');

% Turbine input file
fid = fopen([PathName,FileName, '.fst'], 'wt');
fprintf(fid, '--------------------------------------------------------------------------------\n');
fprintf(fid, '------- FAST INPUT FILE --------------------------------------------------------\n');
fprintf(fid, 'Input files for FAST v7.02.\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- SIMULATION CONTROL --------------------------------------\n');
fprintf(fid, 'False       Echo        - Echo input data to "echo.out" (flag)\n');
fprintf(fid, '   3        ADAMSPrep   - ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create an ADAMS model, 3: do both} (switch)\n');
fprintf(fid, '   1        AnalMode    - Analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch)\n');
fprintf(fid, '   3        NumBl       - Number of blades (-)\n');
fprintf(fid, ' 660.0      TMax        - Total run time (s)\n');
fprintf(fid, '   0.0125   DT          - Integration time step (s)\n');
fprintf(fid, '---------------------- TURBINE CONTROL -----------------------------------------\n');
fprintf(fid, '   0        YCMode      - Yaw control mode {0: none, 1: user-defined from routine UserYawCont, 2: user-defined from Simulink/Labview} (switch)\n');
fprintf(fid, '9999.9      TYCOn       - Time to enable active yaw control (s) [unused when YCMode=0]\n');
fprintf(fid, '   2        PCMode      - Pitch control mode {0: none, 1: user-defined from routine PitchCntrl, 2: user-defined from Simulink/Labview} (switch)\n');
fprintf(fid, '   0.0      TPCOn       - Time to enable active pitch control (s) [unused when PCMode=0]\n');
fprintf(fid, '   3        VSContrl    - Variable-speed control mode {0: none, 1: simple VS, 2: user-defined from routine UserVSCont, 3: user-defined from Simulink/Labview} (switch)\n');
fprintf(fid, ' %5.2f      VS_RtGnSp   - Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]\n', handles.Generator_Speed);
fprintf(fid, ' %5.2f      VS_RtTq     - Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]\n', handles.Generator_Torque);
fprintf(fid, ' %0.8f      VS_Rgn2K    - Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]\n', handles.Generator_K / (60/(2*pi))^2);
fprintf(fid, '  10        VS_SlPc     - Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (percent) [used only when VSContrl=1]\n');
fprintf(fid, '   2        GenModel    - Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]\n');
fprintf(fid, 'True        GenTiStr    - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)\n');
fprintf(fid, 'True        GenTiStp    - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)\n');
fprintf(fid, '9999.9      SpdGenOn    - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]\n');
fprintf(fid, '   0.0      TimGenOn    - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]\n');
fprintf(fid, '9999.9      TimGenOf    - Time to turn off the generator (s) [used only when GenTiStp=True]\n');
fprintf(fid, '   1        HSSBrMode   - HSS brake model {1: simple, 2: user-defined from routine UserHSSBr, 3: user-defined from Labview} (switch)\n');
fprintf(fid, '9999.9      THSSBrDp    - Time to initiate deployment of the HSS brake (s)\n');
fprintf(fid, '9999.9      TiDynBrk    - Time to initiate deployment of the dynamic generator brake [CURRENTLY IGNORED] (s)\n');
fprintf(fid, '9999.9      TTpBrDp(1)  - Time to initiate deployment of tip brake 1 (s)\n');
fprintf(fid, '9999.9      TTpBrDp(2)  - Time to initiate deployment of tip brake 2 (s)\n');
fprintf(fid, '9999.9      TTpBrDp(3)  - Time to initiate deployment of tip brake 3 (s) [unused for 2 blades]\n');
fprintf(fid, '9999.9      TBDepISp(1) - Deployment-initiation speed for the tip brake on blade 1 (rpm)\n');
fprintf(fid, '9999.9      TBDepISp(2) - Deployment-initiation speed for the tip brake on blade 2 (rpm)\n');
fprintf(fid, '9999.9      TBDepISp(3) - Deployment-initiation speed for the tip brake on blade 3 (rpm) [unused for 2 blades]\n');
fprintf(fid, '9999.9      TYawManS    - Time to start override yaw maneuver and end standard yaw control (s)\n');
fprintf(fid, '9999.9      TYawManE    - Time at which override yaw maneuver reaches final yaw angle (s)\n');
fprintf(fid, '   0.0      NacYawF     - Final yaw angle for yaw maneuvers (degrees)\n');
fprintf(fid, '9999.9      TPitManS(1) - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)\n');
fprintf(fid, '9999.9      TPitManS(2) - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)\n');
fprintf(fid, '9999.9      TPitManS(3) - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]\n');
fprintf(fid, '9999.9      TPitManE(1) - Time at which override pitch maneuver for blade 1 reaches final pitch (s)\n');
fprintf(fid, '9999.9      TPitManE(2) - Time at which override pitch maneuver for blade 2 reaches final pitch (s)\n');
fprintf(fid, '9999.9      TPitManE(3) - Time at which override pitch maneuver for blade 3 reaches final pitch (s) [unused for 2 blades]\n');
fprintf(fid, ' %4.3f      BlPitch(1)  - Blade 1 initial pitch (degrees)\n', handles.Blade_PitchOffset);
fprintf(fid, ' %4.3f      BlPitch(2)  - Blade 2 initial pitch (degrees)\n', handles.Blade_PitchOffset);
fprintf(fid, ' %4.3f      BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]\n', handles.Blade_PitchOffset);
fprintf(fid, '   0.0      BlPitchF(1) - Blade 1 final pitch for pitch maneuvers (degrees)\n');
fprintf(fid, '   0.0      BlPitchF(2) - Blade 2 final pitch for pitch maneuvers (degrees)\n');
fprintf(fid, '   0.0      BlPitchF(3) - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]\n');
fprintf(fid, '---------------------- ENVIRONMENTAL CONDITIONS --------------------------------\n');
fprintf(fid, '   9.80665  Gravity     - Gravitational acceleration (m/s^2)\n');
fprintf(fid, '---------------------- FEATURE FLAGS -------------------------------------------\n');
fprintf(fid, 'True        FlapDOF1    - First flapwise blade mode DOF (flag)\n');
fprintf(fid, 'True        FlapDOF2    - Second flapwise blade mode DOF (flag)\n');
fprintf(fid, 'True        EdgeDOF     - First edgewise blade mode DOF (flag)\n');
fprintf(fid, 'False       TeetDOF     - Rotor-teeter DOF (flag) [unused for 3 blades]\n');
fprintf(fid, 'True        DrTrDOF     - Drivetrain rotational-flexibility DOF (flag)\n');
fprintf(fid, 'True        GenDOF      - Generator DOF (flag)\n');
fprintf(fid, 'True        YawDOF      - Yaw DOF (flag)\n');
fprintf(fid, 'True        TwFADOF1    - First fore-aft tower bending-mode DOF (flag)\n');
fprintf(fid, 'True        TwFADOF2    - Second fore-aft tower bending-mode DOF (flag)\n');
fprintf(fid, 'True        TwSSDOF1    - First side-to-side tower bending-mode DOF (flag)\n');
fprintf(fid, 'True        TwSSDOF2    - Second side-to-side tower bending-mode DOF (flag)\n');
fprintf(fid, 'True        CompAero    - Compute aerodynamic forces (flag)\n');
fprintf(fid, 'False       CompNoise   - Compute aerodynamic noise (flag)\n');
fprintf(fid, '---------------------- INITIAL CONDITIONS --------------------------------------\n');
fprintf(fid, '   0.0      OoPDefl     - Initial out-of-plane blade-tip displacement (meters)\n');
fprintf(fid, '   0.0      IPDefl      - Initial in-plane blade-tip deflection (meters)\n');
fprintf(fid, '   0.0      TeetDefl    - Initial or fixed teeter angle (degrees) [unused for 3 blades]\n');
fprintf(fid, '   0.0      Azimuth     - Initial azimuth angle for blade 1 (degrees)\n');
fprintf(fid, '   5.0      RotSpeed    - Initial or fixed rotor speed (rpm)\n');
fprintf(fid, '   0.0      NacYaw      - Initial or fixed nacelle-yaw angle (degrees)\n');
fprintf(fid, '   0.0      TTDspFA     - Initial fore-aft tower-top displacement (meters)\n');
fprintf(fid, '   0.0      TTDspSS     - Initial side-to-side tower-top displacement (meters)\n');
fprintf(fid, '---------------------- TURBINE CONFIGURATION -----------------------------------\n');
fprintf(fid, ' %5.4f      TipRad      - The distance from the rotor apex to the blade tip (meters)\n', handles.Blade_Radius(end));
fprintf(fid, ' %5.4f      HubRad      - The distance from the rotor apex to the blade root (meters)\n', handles.Blade_Radius(1));
fprintf(fid, '   1        PSpnElN     - Number of the innermost blade element which is still part of the pitchable portion of the blade for partial-span pitch control [1 to BldNodes] [CURRENTLY IGNORED] (-)\n');
fprintf(fid, '   0.0      UndSling    - Undersling length [distance from teeter pin to the rotor apex] (meters) [unused for 3 blades]\n');
fprintf(fid, '   0.0      HubCM       - Distance from rotor apex to hub mass [positive downwind] (meters)\n');
fprintf(fid, ' %5.5f      OverHang    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)\n', -handles.Hub_Overhang);
fprintf(fid, ' %5.5f      NacCMxn     - Downwind distance from the tower-top to the nacelle CM (meters)\n', handles.Nacelle_Length*0.7-handles.Hub_Overhang);
fprintf(fid, '   0.0      NacCMyn     - Lateral  distance from the tower-top to the nacelle CM (meters)\n');
fprintf(fid, ' %5.2f      NacCMzn     - Vertical distance from the tower-top to the nacelle CM (meters)\n', handles.Nacelle_Diameter*0.35);
fprintf(fid, ' %5.2f      TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)\n', handles.Tower_Height(end));
fprintf(fid, ' %5.5f      Twr2Shft    - Vertical distance from the tower-top to the rotor shaft (meters)\n', handles.Hub_Height-handles.Tower_Height(end));
fprintf(fid, '   0.0      TwrRBHt     - Tower rigid base height (meters)\n');
fprintf(fid, ' %5.1f      ShftTilt    - Rotor shaft tilt angle (degrees)\n', -handles.Shaft_Tilt);
fprintf(fid, '   0.0      Delta3      - Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]\n');
fprintf(fid, ' %5.1f      PreCone(1)  - Blade 1 cone angle (degrees)\n', -handles.Blade_Cone);
fprintf(fid, ' %5.1f      PreCone(2)  - Blade 2 cone angle (degrees)\n', -handles.Blade_Cone);
fprintf(fid, ' %5.1f      PreCone(3)  - Blade 3 cone angle (degrees) [unused for 2 blades]\n', -handles.Blade_Cone);
fprintf(fid, '   0.0      AzimB1Up    - Azimuth value to use for I/O when blade 1 points up (degrees)\n');
fprintf(fid, '---------------------- MASS AND INERTIA ----------------------------------------\n');
fprintf(fid, '   0.0      YawBrMass   - Yaw bearing mass (kg)\n');
fprintf(fid, ' %5.2E      NacMass     - Nacelle mass (kg)\n', handles.Nacelle_Mass);
fprintf(fid, ' %5.2E      HubMass     - Hub mass (kg)\n', handles.Hub_Mass);
fprintf(fid, '   0.0      TipMass(1)  - Tip-brake mass, blade 1 (kg)\n');
fprintf(fid, '   0.0      TipMass(2)  - Tip-brake mass, blade 2 (kg)\n');
fprintf(fid, '   0.0      TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]\n');
fprintf(fid, '2607.89E3   NacYIner    - Nacelle inertia about yaw axis (kg m^2)\n');
fprintf(fid, ' %5.3f      GenIner     - Generator inertia about HSS (kg m^2)\n', handles.HSS_Inertia);
fprintf(fid, ' 115.926E3  HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)\n');
fprintf(fid, '---------------------- DRIVETRAIN ----------------------------------------------\n');
fprintf(fid, ' %5.1f      GBoxEff     - Gearbox efficiency (percent)\n', 100*handles.Gearbox_Efficiency);
fprintf(fid, ' %5.1f      GenEff      - Generator efficiency [ignored by the Thevenin and user-defined generator models] (percent)\n', 100*handles.Generator_Efficiency);
fprintf(fid, ' %5.1f      GBRatio     - Gearbox ratio (-)\n', handles.Gearbox_Ratio);
fprintf(fid, 'False       GBRevers    - Gearbox reversal {T: if rotor and generator rotate in opposite directions} (flag)\n');
fprintf(fid, '  28.1162E3 HSSBrTqF    - Fully deployed HSS-brake torque (N-m)\n');
fprintf(fid, '   0.6      HSSBrDT     - Time for HSS-brake to reach full deployment once initiated (sec) [used only when HSSBrMode=1]\n');
fprintf(fid, '"Dummy"     DynBrkFi    - File containing a mech-gen-torque vs HSS-speed curve for a dynamic brake [CURRENTLY IGNORED] (quoted string)\n');
fprintf(fid, ' 867.637E6  DTTorSpr    - Drivetrain torsional spring (N-m/rad)\n');
fprintf(fid, '   6.215E6  DTTorDmp    - Drivetrain torsional damper (N-m/(rad/s))\n');
fprintf(fid, '---------------------- SIMPLE INDUCTION GENERATOR ------------------------------\n');
fprintf(fid, '9999.9      SIG_SlPc    - Rated generator slip percentage (percent) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '9999.9      SIG_SySp    - Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '9999.9      SIG_RtTq    - Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '9999.9      SIG_PORt    - Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]\n');
fprintf(fid, '---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR -----------------\n');
fprintf(fid, '9999.9      TEC_Freq    - Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9998        TEC_NPol    - Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9      TEC_SRes    - Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9      TEC_RRes    - Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9      TEC_VLL     - Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9      TEC_SLR     - Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9      TEC_RLR     - Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '9999.9      TEC_MR      - Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]\n');
fprintf(fid, '---------------------- PLATFORM ------------------------------------------------\n');
fprintf(fid, '   0        PtfmModel   - Platform model {0: none, 1: onshore, 2: fixed bottom offshore, 3: floating offshore} (switch)\n');
fprintf(fid, '"Dummy"     PtfmFile    - Name of file containing platform properties (quoted string) [unused when PtfmModel=0]\n');
fprintf(fid, '---------------------- TOWER ---------------------------------------------------\n');
fprintf(fid, '  %i        TwrNodes    - Number of tower nodes used for analysis (-)\n', length(handles.Tower_Height));
fprintf(fid, '"%s"          TwrFile     - Name of file containing tower properties (quoted string)\n', [FileName, '_tower.dat']);
fprintf(fid, '---------------------- NACELLE-YAW ---------------------------------------------\n');
fprintf(fid, '9028.32E6   YawSpr      - Nacelle-yaw spring constant (N-m/rad)\n');
fprintf(fid, '  19.16E6   YawDamp     - Nacelle-yaw damping constant (N-m/(rad/s))\n');
fprintf(fid, '   0.0      YawNeut     - Neutral yaw position--yaw spring force is zero at this yaw (degrees)\n');
fprintf(fid, '---------------------- FURLING -------------------------------------------------\n');
fprintf(fid, 'False       Furling     - Read in additional model properties for furling turbine (flag)\n');
fprintf(fid, '"Dummy"     FurlFile    - Name of file containing furling properties (quoted string) [unused when Furling=False]\n');
fprintf(fid, '---------------------- ROTOR-TEETER --------------------------------------------\n');
fprintf(fid, '   0        TeetMod     - Rotor-teeter spring/damper model {0: none, 1: standard, 2: user-defined from routine UserTeet} (switch) [unused for 3 blades]\n');
fprintf(fid, '   0.0      TeetDmpP    - Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '   0.0      TeetDmp     - Rotor-teeter damping constant (N-m/(rad/s)) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '   0.0      TeetCDmp    - Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '   0.0      TeetSStP    - Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '   0.0      TeetHStP    - Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '   0.0      TeetSSSp    - Rotor-teeter soft-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '   0.0      TeetHSSp    - Rotor-teeter hard-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]\n');
fprintf(fid, '---------------------- TIP-BRAKE -----------------------------------------------\n');
fprintf(fid, '   0.0      TBDrConN    - Tip-brake drag constant during normal operation, Cd*Area (m^2)\n');
fprintf(fid, '   0.0      TBDrConD    - Tip-brake drag constant during fully-deployed operation, Cd*Area (m^2)\n');
fprintf(fid, '   0.0      TpBrDT      - Time for tip-brake to reach full deployment once released (sec)\n');
fprintf(fid, '---------------------- BLADE ---------------------------------------------------\n');
fprintf(fid, '"%s"                  BldFile(1)  - Name of file containing properties for blade 1 (quoted string)\n', [FileName, '_blade.dat']);
fprintf(fid, '"%s"                  BldFile(2)  - Name of file containing properties for blade 2 (quoted string)\n', [FileName, '_blade.dat']);
fprintf(fid, '"%s"                  BldFile(3)  - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]\n', [FileName, '_blade.dat']);
fprintf(fid, '---------------------- AERODYN -------------------------------------------------\n');
fprintf(fid, '"%s"                ADFile      - Name of file containing AeroDyn input parameters (quoted string)\n', [FileName, '_AeroDyn.ipt']);
fprintf(fid, '---------------------- NOISE ---------------------------------------------------\n');
fprintf(fid, '"Dummy"     NoiseFile   - Name of file containing aerodynamic noise input parameters (quoted string) [used only when CompNoise=True]\n');
fprintf(fid, '---------------------- ADAMS ---------------------------------------------------\n');
fprintf(fid, '"%s"          ADAMSFile   - Name of file containing ADAMS-specific input parameters (quoted string) [unused when ADAMSPrep=1]\n', [FileName, '_ADAMS.dat']);
fprintf(fid, '---------------------- LINEARIZATION CONTROL -----------------------------------\n');
fprintf(fid, '"Dummy"                 LinFile     - Name of file containing FAST linearization parameters (quoted string) [unused when AnalMode=1]\n');
fprintf(fid, '---------------------- OUTPUT --------------------------------------------------\n');
fprintf(fid, 'True        SumPrint    - Print summary data to "<RootName>.fsm" (flag)\n');
fprintf(fid, '1           OutFileFmt  - Format for tabular (time-marching) output file(s) (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both) (switch)\n');
fprintf(fid, 'True        TabDelim    - Use tab delimiters in text tabular output file? (flag)\n');
fprintf(fid, '"ES10.3E2"  OutFmt      - Format used for text tabular output (except time).  Resulting field should be 10 characters. (quoted string)  [not checked for validity!]\n');
fprintf(fid, '   0.0      TStart      - Time to begin tabular output (s)\n');
fprintf(fid, '   1        DecFact     - Decimation factor for tabular output {1: output every time step} (-)\n');
fprintf(fid, '   1.0      SttsTime    - Amount of time between screen status messages (sec)\n');
fprintf(fid, '  -3.09528  NcIMUxn     - Downwind distance from the tower-top to the nacelle IMU (meters)\n');
fprintf(fid, '   0.0      NcIMUyn     - Lateral  distance from the tower-top to the nacelle IMU (meters)\n');
fprintf(fid, '   2.23336  NcIMUzn     - Vertical distance from the tower-top to the nacelle IMU (meters)\n');
fprintf(fid, '   1.912    ShftGagL    - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)\n');
fprintf(fid, '   0        NTwGages    - Number of tower nodes that have strain gages for output [0 to 9] (-)\n');
fprintf(fid, '            TwrGagNd    - List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]\n');
fprintf(fid, '   3        NBlGages    - Number of blade nodes that have strain gages for output [0 to 9] (-)\n');
fprintf(fid, '5,9,13     BldGagNd    - List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]\n');
fprintf(fid, '            OutList     - The next line(s) contains a list of output parameters.  See OutList.txt for a listing of available output channels, (-)\n');
fprintf(fid, '"uWind	, vWind,   , wWind"\n');
fprintf(fid, '"OoPDefl1 , OoPDefl2 , OoPDefl3"\n');
fprintf(fid, '"IPDefl1  , IPDefl2  , IPDefl3"\n');
fprintf(fid, '"NcIMUTAxs, NcIMUTAys, NcIMUTAzs"\n');
fprintf(fid, '"RootMOoP1, RootMOoP2, RootMOoP3"\n');
fprintf(fid, '"RootMIP1 , RootMIP2 , RootMIP3"\n');
fprintf(fid, '"RootMFlp1, RootMFlp2, RootMFlp3"\n');
fprintf(fid, '"RootMEdg1, RootMEdg2, RootMEdg3"\n');
fprintf(fid, '"TwrBsMxt , TwrBsMyt , TwrBsMzt"\n');
fprintf(fid, '"LSSTipMya, LSSTipMza, LSSTipVxa"\n');
fprintf(fid, '"GenPwr   , GenTq    , GenSpeed"\n');
fprintf(fid, '"BlPitch1 , Azimuth"\n');
fprintf(fid, '"RotPwr   , RotThrust"\n');
fprintf(fid, '"RotSpeed , HSShftTq"\n');
fprintf(fid, 'END of FAST input file (the word "END" must appear in the first 3 columns of this last line).\n');
fprintf(fid, '--------------------------------------------------------------------------------\n');
fclose(fid);

% Tower input file
fid = fopen([PathName,FileName, '_tower.dat'], 'wt');
fprintf(fid, '--------------------------------------------------------------------------------\n');
fprintf(fid, '---------------------- FAST TOWER FILE -----------------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- TOWER PARAMETERS ----------------------------------------\n');
fprintf(fid, '  %i        NTwInpSt    - Number of input stations to specify tower geometry\n', length(handles.Tower_Height));
fprintf(fid, 'False       CalcTMode   - Calculate tower mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)\n');
fprintf(fid, '   1.0      TwrFADmp(1) - Tower 1st fore-aft mode structural damping ratio (percent)\n');
fprintf(fid, '   1.0      TwrFADmp(2) - Tower 2nd fore-aft mode structural damping ratio (percent)\n');
fprintf(fid, '   1.0      TwrSSDmp(1) - Tower 1st side-to-side mode structural damping ratio (percent)\n');
fprintf(fid, '   1.0      TwrSSDmp(2) - Tower 2nd side-to-side mode structural damping ratio (percent)\n');
fprintf(fid, '---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------\n');
fprintf(fid, '   1.0      FAStTunr(1) - Tower fore-aft modal stiffness tuner, 1st mode (-)\n');
fprintf(fid, '   1.0      FAStTunr(2) - Tower fore-aft modal stiffness tuner, 2nd mode (-)\n');
fprintf(fid, '   1.0      SSStTunr(1) - Tower side-to-side stiffness tuner, 1st mode (-)\n');
fprintf(fid, '   1.0      SSStTunr(2) - Tower side-to-side stiffness tuner, 2nd mode (-)\n');
fprintf(fid, '   1.0      AdjTwMa     - Factor to adjust tower mass density (-)\n');
fprintf(fid, '   1.0      AdjFASt     - Factor to adjust tower fore-aft stiffness (-)\n');
fprintf(fid, '   1.0      AdjSSSt     - Factor to adjust tower side-to-side stiffness (-)\n');
fprintf(fid, '---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------\n');
fprintf(fid, 'HtFract  TMassDen  TwFAStif   TwSSStif   TwGJStif   TwEAStif   TwFAIner  TwSSIner  TwFAcgOf  TwSScgOf\n');
fprintf(fid, '(-)      (kg/m)    (Nm^2)     (Nm^2)     (Nm^2)     (N)        (kg m)    (kg m)    (m)       (m)\n');
for i = 1:length(handles.Tower_Height)
    fprintf(fid, '%0.2f\t%7.2f\t%7.3E\t%7.3E\t%7.3E\t%7.3E\t%7.1f\t%7.1f\t%7.1f\t%7.1f\n', ...
        handles.Tower_Height(i)/handles.Tower_Height(end), ...
        handles.Tower_Mass(i), ...
        handles.Tower_EI(i), ...
        handles.Tower_EI(i), ...
        handles.Tower_GJ(i), ...
        handles.Tower_EA(i), ...
        handles.Tower_Iner(i), ...
        handles.Tower_Iner(i), ...
        0, ...
        0);
end
fprintf(fid, '---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------\n');
fprintf(fid, '    %9.4f   TwFAM1Sh(2) - Mode 1, coefficient of x^2 term\n', handles.Tower_ForeAft1_coeff(5));
fprintf(fid, '    %9.4f   TwFAM1Sh(3) -       , coefficient of x^3 term\n', handles.Tower_ForeAft1_coeff(4));
fprintf(fid, '    %9.4f   TwFAM1Sh(4) -       , coefficient of x^4 term\n', handles.Tower_ForeAft1_coeff(3));
fprintf(fid, '    %9.4f   TwFAM1Sh(5) -       , coefficient of x^5 term\n', handles.Tower_ForeAft1_coeff(2));
fprintf(fid, '    %9.4f   TwFAM1Sh(6) -       , coefficient of x^6 term\n', handles.Tower_ForeAft1_coeff(1));
fprintf(fid, '    %9.4f   TwFAM2Sh(2) - Mode 2, coefficient of x^2 term\n', handles.Tower_ForeAft2_coeff(5));
fprintf(fid, '    %9.4f   TwFAM2Sh(3) -       , coefficient of x^3 term\n', handles.Tower_ForeAft2_coeff(4));
fprintf(fid, '    %9.4f   TwFAM2Sh(4) -       , coefficient of x^4 term\n', handles.Tower_ForeAft2_coeff(3));
fprintf(fid, '    %9.4f   TwFAM2Sh(5) -       , coefficient of x^5 term\n', handles.Tower_ForeAft2_coeff(2));
fprintf(fid, '    %9.4f   TwFAM2Sh(6) -       , coefficient of x^6 term\n', handles.Tower_ForeAft2_coeff(1));
fprintf(fid, '---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------\n');
fprintf(fid, '    %9.4f   TwSSM1Sh(2) - Mode 1, coefficient of x^2 term\n', handles.Tower_SideSide1_coeff(5));
fprintf(fid, '    %9.4f   TwSSM1Sh(3) -       , coefficient of x^3 term\n', handles.Tower_SideSide1_coeff(4));
fprintf(fid, '    %9.4f   TwSSM1Sh(4) -       , coefficient of x^4 term\n', handles.Tower_SideSide1_coeff(3));
fprintf(fid, '    %9.4f   TwSSM1Sh(5) -       , coefficient of x^5 term\n', handles.Tower_SideSide1_coeff(2));
fprintf(fid, '    %9.4f   TwSSM1Sh(6) -       , coefficient of x^6 term\n', handles.Tower_SideSide1_coeff(1));
fprintf(fid, '    %9.4f   TwSSM2Sh(2) - Mode 2, coefficient of x^2 term\n', handles.Tower_SideSide2_coeff(5));
fprintf(fid, '    %9.4f   TwSSM2Sh(3) -       , coefficient of x^3 term\n', handles.Tower_SideSide2_coeff(4));
fprintf(fid, '    %9.4f   TwSSM2Sh(4) -       , coefficient of x^4 term\n', handles.Tower_SideSide2_coeff(3));
fprintf(fid, '    %9.4f   TwSSM2Sh(5) -       , coefficient of x^5 term\n', handles.Tower_SideSide2_coeff(2));
fprintf(fid, '    %9.4f   TwSSM2Sh(6) -       , coefficient of x^6 term\n', handles.Tower_SideSide2_coeff(1));
fclose(fid);

% Blade input file
fid = fopen([PathName,FileName, '_blade.dat'], 'wt');
fprintf(fid, '--------------------------------------------------------------------------------\n');
fprintf(fid, '---------------------- FAST INDIVIDUAL BLADE FILE ------------------------------\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- BLADE PARAMETERS ----------------------------------------\n');
fprintf(fid, '  %i        NBlInpSt    - Number of blade input stations (-)\n', length(handles.Blade_Radius));
fprintf(fid, 'False       CalcBMode   - Calculate blade mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)\n');
fprintf(fid, '   0.477465 BldFlDmp(1) - Blade flap mode #1 structural damping in percent of critical (percent)\n');
fprintf(fid, '   0.477465 BldFlDmp(2) - Blade flap mode #2 structural damping in percent of critical (percent)\n');
fprintf(fid, '   0.477465 BldEdDmp(1) - Blade edge mode #1 structural damping in percent of critical (percent)\n');
fprintf(fid, '---------------------- BLADE ADJUSTMENT FACTORS --------------------------------\n');
fprintf(fid, '   1.0      FlStTunr(1) - Blade flapwise modal stiffness tuner, 1st mode (-)\n');
fprintf(fid, '   1.0      FlStTunr(2) - Blade flapwise modal stiffness tuner, 2nd mode (-)\n');
fprintf(fid, '   1.04536  AdjBlMs     - Factor to adjust blade mass density (-)\n');
fprintf(fid, '   1.0      AdjFlSt     - Factor to adjust blade flap stiffness (-)\n');
fprintf(fid, '   1.0      AdjEdSt     - Factor to adjust blade edge stiffness (-)\n');
fprintf(fid, '---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------\n');
fprintf(fid, 'BlFract  AeroCent  StrcTwst  BMassDen  FlpStff     EdgStff     GJStff     EAStff      Alpha  FlpIner  EdgIner  PrecrvRef  PreswpRef  FlpcgOf  EdgcgOf   FlpEAOf  EdgEAOf\n');
fprintf(fid, '(-)      (-)       (deg)     (kg/m)    (Nm^2)      (Nm^2)      (Nm^2)     (N)         (-)    (kg m)   (kg m)   (m)        (m)        (m)      (m)       (m)      (m)\n');
for i = 1:length(handles.Blade_Radius)
    fprintf(fid, '%7.5f\t%7.5f\t%7.3f\t%7.3f\t%7.5E\t%7.5E\t%7.5E\t%7.5E\t%7.1f\t%7.2f\t%7.2f\t%7.1f\t%7.1f\t%7.1f\t%7.5f\t%7.1f\t%7.1f\n', ...
        (handles.Blade_Radius(i)-handles.Blade_Radius(1))/(handles.Blade_Radius(end)-handles.Blade_Radius(1)), ...
        handles.Blade_ac(i), ...
        handles.Blade_Twist(i), ...
        handles.Blade_Mass(i), ...
        handles.Blade_EIflap(i), ...
        handles.Blade_EIedge(i), ...
        handles.Blade_GJ(i), ...
        handles.Blade_EA(i), ...
        0, ...
        handles.Blade_FlapIner(i), ...
        handles.Blade_EdgeIner(i), ...
        0, ...
        0, ...
        0, ...
        handles.Blade_cg(i), ...
        0, ...
        0);
end
fprintf(fid, '---------------------- BLADE MODE SHAPES ---------------------------------------\n');
fprintf(fid, '    %9.4f   BldFl1Sh(2) - Flap mode 1, coeff of x^2\n', handles.Blade_Flap1_coeff(5));
fprintf(fid, '    %9.4f   BldFl1Sh(3) -            , coeff of x^3\n', handles.Blade_Flap1_coeff(4));
fprintf(fid, '    %9.4f   BldFl1Sh(4) -            , coeff of x^4\n', handles.Blade_Flap1_coeff(3));
fprintf(fid, '    %9.4f   BldFl1Sh(5) -            , coeff of x^5\n', handles.Blade_Flap1_coeff(2));
fprintf(fid, '    %9.4f   BldFl1Sh(6) -            , coeff of x^6\n', handles.Blade_Flap1_coeff(1));
fprintf(fid, '    %9.4f   BldFl2Sh(2) - Flap mode 2, coeff of x^2\n', handles.Blade_Flap2_coeff(5));
fprintf(fid, '    %9.4f   BldFl2Sh(3) -            , coeff of x^3\n', handles.Blade_Flap2_coeff(4));
fprintf(fid, '    %9.4f   BldFl2Sh(4) -            , coeff of x^4\n', handles.Blade_Flap2_coeff(3));
fprintf(fid, '    %9.4f   BldFl2Sh(5) -            , coeff of x^5\n', handles.Blade_Flap2_coeff(2));
fprintf(fid, '    %9.4f   BldFl2Sh(6) -            , coeff of x^6\n', handles.Blade_Flap2_coeff(1));
fprintf(fid, '    %9.4f   BldEdgSh(2) - Edge mode 1, coeff of x^2\n', handles.Blade_Edge1_coeff(5));
fprintf(fid, '    %9.4f   BldEdgSh(3) -            , coeff of x^3\n', handles.Blade_Edge1_coeff(4));
fprintf(fid, '    %9.4f   BldEdgSh(4) -            , coeff of x^4\n', handles.Blade_Edge1_coeff(3));
fprintf(fid, '    %9.4f   BldEdgSh(5) -            , coeff of x^5\n', handles.Blade_Edge1_coeff(2));
fprintf(fid, '    %9.4f   BldEdgSh(6) -            , coeff of x^6\n', handles.Blade_Edge1_coeff(1));
fclose(fid);

% Airfoil files
copyfile([pwd, '\subfunctions\AeroData\Cylinder1.dat'], [PathName, FileName, '_airfoil_Cylinder1.dat']);
copyfile([pwd, '\subfunctions\AeroData\Cylinder2.dat'], [PathName, FileName, '_airfoil_Cylinder2.dat']);
copyfile([pwd, '\subfunctions\AeroData\DU40_A17.dat'], [PathName, FileName, '_airfoil_DU40_A17.dat']);
copyfile([pwd, '\subfunctions\AeroData\DU35_A17.dat'], [PathName, FileName, '_airfoil_DU35_A17.dat']);
copyfile([pwd, '\subfunctions\AeroData\DU30_A17.dat'], [PathName, FileName, '_airfoil_DU30_A17.dat']);
copyfile([pwd, '\subfunctions\AeroData\DU25_A17.dat'], [PathName, FileName, '_airfoil_DU25_A17.dat']);
copyfile([pwd, '\subfunctions\AeroData\DU21_A17.dat'], [PathName, FileName, '_airfoil_DU21_A17.dat']);
copyfile([pwd, '\subfunctions\AeroData\NACA64_A17.dat'], [PathName, FileName, '_airfoil_NACA64_A17.dat']);

% AeroDyn input file
fid = fopen([PathName,FileName, '_AeroDyn.ipt'], 'wt');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, 'SI          SysUnits    - System of units for used for input and output [must be SI for FAST] (unquoted string)\n');
fprintf(fid, 'BEDDOES     StallMod    - Dynamic stall included [BEDDOES or STEADY] (unquoted string)\n');
fprintf(fid, 'USE_CM      UseCm       - Use aerodynamic pitching moment model? [USE_CM or NO_CM] (unquoted string)\n');
fprintf(fid, 'EQUIL       InfModel    - Inflow model [DYNIN or EQUIL] (unquoted string)\n');
fprintf(fid, 'SWIRL       IndModel    - Induction-factor model [NONE or WAKE or SWIRL] (unquoted string)\n');
fprintf(fid, '   0.005    AToler      - Induction-factor tolerance (convergence criteria) (-)\n');
fprintf(fid, 'PRANDtl     TLModel     - Tip-loss model (EQUIL only) [PRANDtl, GTECH, or NONE] (unquoted string)\n');
fprintf(fid, 'PRANDtl     HLModel     - Hub-loss model (EQUIL only) [PRANdtl or NONE] (unquoted string)\n');
fprintf(fid, '"wind.wnd"    WindFile    - Name of file containing wind data (quoted string)\n');
fprintf(fid, ' %5.1f      HH          - Wind reference (hub) height [TowerHt+Twr2Shft+OverHang*SIN(ShftTilt)] (m)\n', handles.Hub_Height);
fprintf(fid, '   0.0      TwrShad     - Tower-shadow velocity deficit (-)\n');
fprintf(fid, '9999.9      ShadHWid    - Tower-shadow half width (m)\n');
fprintf(fid, '9999.9      T_Shad_Refpt - Tower-shadow reference point (m)\n');
fprintf(fid, '   1.225    AirDens     - Air density (kg/m^3)\n');
fprintf(fid, '   1.464E-5 KinVisc     - Kinematic air viscosity [CURRENTLY IGNORED] (m^2/sec)\n');
fprintf(fid, '   0.02479  DTAero      - Time interval for aerodynamic calculations (sec)\n');
fprintf(fid, '   8        NumFoil     - Number of airfoil files (-)\n');
fprintf(fid, '"%s"                         FoilNm      - Names of the airfoil files [NumFoil lines] (quoted strings)\n', [FileName, '_airfoil_Cylinder1.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_Cylinder2.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_DU40_A17.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_DU35_A17.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_DU30_A17.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_DU25_A17.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_DU21_A17.dat']);
fprintf(fid, '"%s"\n', [FileName, '_airfoil_NACA64_A17.dat']);
fprintf(fid, '  %i        BldNodes    - Number of blade nodes used for analysis (-)\n', length(Blade_r));
fprintf(fid, 'RNodes   AeroTwst  DRNodes  Chord  NFoil  PrnElm\n');
for i = 1:length(Blade_r)
    fprintf(fid, '%7.4f\t%7.3f\t%7.4f\t%7.3f\t%i\tNOPRINT\n', ...
        Blade_r(i), ...
        Blade_Twist(i), ...
        Blade_dr(i), ...
        Blade_Chord(i), ...
        Blade_NFoil(i));
end
fclose(fid);

% ADAMS input file
fid = fopen([PathName,FileName, '_ADAMS.dat'], 'wt');
fprintf(fid, '--------------------------------------------------------------------------------\n');
fprintf(fid, '---------------------- FAST 2 ADAMS PREPROCESSOR, ADAMS-SPECIFIC DATA FILE -----\n');
fprintf(fid, 'Created %s.\n', datestr(now));
fprintf(fid, '---------------------- FEATURE FLAGS -------------------------------------------\n');
fprintf(fid, 'True        SaveGrphcs  - Save GRAPHICS output (flag)\n');
fprintf(fid, 'True        MakeLINacf  - Make an ADAMS/LINEAR control / command file (flag)\n');
fprintf(fid, '---------------------- DAMPING PARAMETERS --------------------------------------\n');
fprintf(fid, '   0.01     CRatioTGJ   - Ratio of damping to stiffness for the tower torsion     deflection  (-)\n');
fprintf(fid, '   0.01     CRatioTEA   - Ratio of damping to stiffness for the tower extensional deflection  (-)\n');
fprintf(fid, '   0.01     CRatioBGJ   - Ratio of damping to stiffness for the blade torsion     deflections (-)\n');
fprintf(fid, '   0.01     CRatioBEA   - Ratio of damping to stiffness for the blade extensional deflections (-)\n');
fprintf(fid, '---------------------- BLADE PITCH ACTUATOR PARAMETERS -------------------------\n');
fprintf(fid, ' 971.350E6  BPActrSpr   - Blade pitch actuator spring stiffness constant (N-m/rad)\n');
fprintf(fid, '   0.206E6  BPActrDmp   - Blade pitch actuator damping          constant (N-m/(rad/s))\n');
fprintf(fid, '---------------------- GRAPHICS PARAMETERS -------------------------------------\n');
fprintf(fid, '  20        NSides      - Number of sides used in GRAPHICS CYLINDER and FRUSTUM statements (-)\n');
fprintf(fid, '   3.000    TwrBaseRad  - Tower base radius used for linearly tapered tower GRAPHICS CYLINDERs (m)\n');
fprintf(fid, '   1.935    TwrTopRad   - Tower top  radius used for linearly tapered tower GRAPHICS CYLINDERs (m)\n');
fprintf(fid, '   7.0      NacLength   - Length of nacelle used for the nacelle GRAPHICS (m)\n');
fprintf(fid, '   1.75     NacRadBot   - Bottom (opposite rotor) radius of nacelle FRUSTUM used for the nacelle GRAPHICS (m)\n');
fprintf(fid, '   1.75     NacRadTop   - Top    (rotor end)      radius of nacelle FRUSTUM used for the nacelle GRAPHICS (m)\n');
fprintf(fid, '   1.0      GBoxLength  - Length, width, and height of the gearbox BOX for gearbox GRAPHICS (m)\n');
fprintf(fid, '   2.39     GenLength   - Length of the generator CYLINDER used for generator GRAPHICS (m)\n');
fprintf(fid, '   1.195    HSSLength   - Length of the high-speed shaft CYLINDER used for HSS GRAPHICS (m)\n');
fprintf(fid, '   4.78     LSSLength   - Length of the low-speed shaft CYLINDER used for LSS GRAPHICS (m)\n');
fprintf(fid, '   0.75     GenRad      - Radius of the generator CYLINDER used for generator GRAPHICS (m)\n');
fprintf(fid, '   0.2      HSSRad      - Radius of the high-speed shaft CYLINDER used for HSS GRAPHICS (m)\n');
fprintf(fid, '   0.4      LSSRad      - Radius of the low -speed shaft CYLINDER used for LSS GRAPHICS (m)\n');
fprintf(fid, '   0.875    HubCylRad   - Radius of hub CYLINDER used for hub GRAPHICS (m)\n');
fprintf(fid, '   0.18     ThkOvrChrd  - Ratio of blade thickness to blade chord used for blade element BOX GRAPHICS (-)\n');
fprintf(fid, '   0.0      BoomRad     - Radius of the tail boom CYLINDER used for tail boom GRAPHICS (m)\n');
fclose(fid);

% Enable buttons
set(handles.Linearize, 'Enable', 'on')
set(handles.WindSpeed_From, 'Enable', 'on')
set(handles.WindSpeed_To, 'Enable', 'on')
set(handles.WindSpeed_Step, 'Enable', 'on')
set(handles.WindSpeed_From, 'String', '5')
set(handles.WindSpeed_To, 'String', '25')
set(handles.WindSpeed_Step, 'String', '1')

%% Linearization
function Linearize_Callback(hObject, eventdata, handles)

% Wind speed range
WindSpeeds = str2double(get(handles.WindSpeed_From, 'String')):str2double(get(handles.WindSpeed_Step, 'String')):str2double(get(handles.WindSpeed_To, 'String'));
if sum(WindSpeeds) == 0 || isnan(sum(WindSpeeds))
    errordlg('Invalid wind speed range.', 'Error')
else
    
set(hObject, 'Enable', 'off');

disp(' Starting linearization...')
disp(' ')

% Get geometry from handles
Blade_Radius = handles.Blade_Radius;
Blade_Twist = handles.Blade_Twist * pi/180;
Blade_Chord = handles.Blade_Chord;
Blade_NFoil = handles.Blade_NFoil;
NBlades = handles.NBlades;
WindSpeed_Cutin = handles.WindSpeed_Cutin;
WindSpeed_Cutout = handles.WindSpeed_Cutout;
Generator_Efficiency = handles.Generator_Efficiency;
Generator_Speed = handles.Generator_Speed;
Gearbox_Efficiency = handles.Gearbox_Efficiency;
Gearbox_Ratio = handles.Gearbox_Ratio;
Rated_Power = handles.Rated_Power;
TSR = handles.Rated_TipSpeedRatio;
AoA = handles.AoA * pi/180;
Pitch = handles.Blade_PitchOffset * pi/180;
for i = 1:length(Blade_NFoil)
    CL(:,i) = handles.CL(:,Blade_NFoil(i));
    CD(:,i) = handles.CD(:,Blade_NFoil(i));
end

% Find maximum CP and optimum TSR
TSR = 0:0.1:20;
CP = zeros(size(TSR));
for i = 2:length(TSR)

    % Local tip speed ratio and solidity
    r = Blade_Radius/Blade_Radius(end);
    lambdar = TSR(i) * r;
    sigmar = NBlades*Blade_Chord./(2*pi*Blade_Radius);

    % Initial induction factors
    anew = 1/3*ones(size(Blade_Radius));
    a_new = zeros(size(Blade_Radius));
    for iter = 1:100

        a = anew;
        a_ = a_new;

        % Inflow angle
        phi = real(atan((1-a)./((1+a_).*lambdar)));

        % Tip loss correction
        F = 2/pi * acos(exp(-NBlades/2*(1-r)./(r.*sin(phi))));

        % Aerodynamic force coefficients
        alpha = phi - Pitch - Blade_Twist;
        for j = 1:length(Blade_Radius)
            Cl(j) = interp1(AoA,CL(:,j),alpha(j));
            Cd(j) = interp1(AoA,CD(:,j),alpha(j));
        end
        Cl = Cl(:);
        Cd = Cd(:);
        Cl(isnan(Cl)) = 1e-6;
        Cd(isnan(Cd)) = 0;
        Cl(Cl == 0) = 1e-6;

        % Thrust coefficient
        CT = sigmar.*(1-a).^2.*(Cl.*cos(phi)+Cd.*sin(phi))./(sin(phi).^2);

        % New induction factors
        for j = 1:length(Blade_Radius)
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

% Estimate for rated wind speed
Urated = (Rated_Power/(0.5 * 1.225 * pi*Blade_Radius(end)^2 * max(CP) * Gearbox_Efficiency * Generator_Efficiency))^(1/3);

% Rotor speeds
TSR = TSR(CP==max(CP));
TSR = TSR(1);
RPM = TSR*WindSpeeds/Blade_Radius(end) * 60/(2*pi);
RPM(RPM > (Generator_Speed/Gearbox_Ratio)) = Generator_Speed / Gearbox_Ratio;

% Avoiding some common errors
% (1) Set the gearbox efficiency (to avoid error that ADAMS cannot handle
% nonideal gearboxes)
Gearbox_Efficiency = 1;

% (2) Scale the HSS inertia from the reference machine (to avoid issues
% with reaching very high rotor speeds during the linearization of some
% direct-drive machines). Assuming that P ~ Mass, Radius² ~ P^(2/3), we get
% a relation in the shape of I/Iref = (P/Pref)^(5/3) * (RPM_ref/RPM)^2.
HSS_Inertia = 534.116 * (Rated_Power/(5e6))^(5/3) * (12.1*97/Generator_Speed)^2;

% Copy other input files
copyfile(get(handles.TowerFile_textbox, 'String'), [pwd, '\subfunctions\linearization_tower.dat']);
copyfile(get(handles.BladeFile_textbox, 'String'), [pwd, '\subfunctions\linearization_blade.dat']);
copyfile(get(handles.ADAMSFile_textbox, 'String'), [pwd, '\subfunctions\linearization_ADAMS.dat']);

% Aerodyn input file
fid = fopen(get(handles.AeroDynFile_textbox, 'String'),'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
A{3} = 'STEADY     StallMod    - Dynamic stall included [BEDDOES or STEADY] (unquoted string)';
A{7} = '   1E-6    AToler      - Induction-factor tolerance (convergence criteria) (-)';
A{10} = '"steadywind.hh"    WindFile    - Name of file containing wind data (quoted string)\n';
A{19} = '"AeroData\Cylinder1.dat"                         FoilNm      - Names of the airfoil files [NumFoil lines] (quoted strings)';
A{20} = '"AeroData\Cylinder2.dat"';
A{21} = '"AeroData\DU40_A17.dat"';
A{22} = '"AeroData\DU35_A17.dat"';
A{23} = '"AeroData\DU30_A17.dat"';
A{24} = '"AeroData\DU25_A17.dat"';
A{25} = '"AeroData\DU21_A17.dat"';
A{26} = '"AeroData\NACA64_A17.dat"';
fid = fopen([pwd, '\subfunctions\linearization_AeroDyn.ipt'], 'wt');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);

% Run linearization
sysm = cell(length(WindSpeeds),1);
for j = 1:length(WindSpeeds)
    
    % Select trim case
    if WindSpeeds(j) <= Urated
        TrimCase = 2;
    else
        TrimCase = 3;
    end
    
    % Turbine input file
    fid = fopen(get(handles.TurbineFile_textbox, 'String'),'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    A{8} = '   2        AnalMode    - Analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch)';
    A{10} = '1500.0      TMax        - Total run time (s)';
    A{15} = '   0        PCMode      - Pitch control mode {0: none, 1: user-defined from routine PitchCntrl, 2: user-defined from Simulink/Labview} (switch)';
    A{17} = '   1        VSContrl    - Variable-speed control mode {0: none, 1: simple VS, 2: user-defined from routine UserVSCont, 3: user-defined from Simulink/Labview} (switch)';
    A{73} = ['  ', num2str(RPM(j), '%2.2f'), '      RotSpeed    - Initial or fixed rotor speed (rpm)'];
    A{104} = [' ', num2str(HSS_Inertia, '%5.3f'), '      GenIner     - Generator inertia about HSS (kg m^2)'];
    A{107} = [' ', num2str(100*Gearbox_Efficiency, '%5.1f'), '      GBoxEff     - Gearbox efficiency (percent)'];
    A{135} = ['"', pwd, '\subfunctions\linearization_tower.dat"          TwrFile     - Name of file containing tower properties (quoted string)'];
    A{157} = ['"', pwd, '\subfunctions\linearization_blade.dat"                  BldFile(1)  - Name of file containing properties for blade 1 (quoted string)'];
    A{158} = ['"', pwd, '\subfunctions\linearization_blade.dat"                  BldFile(2)  - Name of file containing properties for blade 2 (quoted string)'];
    A{159} = ['"', pwd, '\subfunctions\linearization_blade.dat"                  BldFile(3)  - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]'];
    A{161} = ['"', pwd, '\subfunctions\linearization_Aerodyn.ipt"                 LinFile     - Name of file containing FAST linearization parameters (quoted string) [unused when AnalMode=1]'];
    A{165} = ['"', pwd, '\subfunctions\linearization_ADAMS.dat"          ADAMSFile   - Name of file containing ADAMS-specific input parameters (quoted string) [unused when ADAMSPrep=1]'];
    A{167} = ['"', pwd, '\subfunctions\linearization_linear.dat"                 LinFile     - Name of file containing FAST linearization parameters (quoted string) [unused when AnalMode=1]'];
    fid = fopen([pwd, '\subfunctions\linearization.fst'], 'wt');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);
    
    % Linearization input file
    fid = fopen([pwd, '\subfunctions\linearization_linear.dat'], 'wt');
    fprintf(fid, '--------------------------------------------------------------------------------\n');
    fprintf(fid, '---------------------- FAST LINEARIZATION CONTROL FILE -------------------------\n');
    fprintf(fid, 'Created %s.\n', datestr(now));
    fprintf(fid, '---------------------- PERIODIC STEADY STATE SOLUTION --------------------------\n');
    fprintf(fid, 'True        CalcStdy    - Calculate periodic steady state condition {False: linearize about initial conditions} (flag)\n');
    fprintf(fid, '   %i        TrimCase    - Trim case {1: find nacelle yaw, 2: find generator torque, 3: find collective blade pitch} (switch) [used only when CalcStdy=True and GenDOF=True]\n', TrimCase);
    fprintf(fid, '   0.0001   DispTol     - Convergence tolerance for the 2-norm of displacements in the periodic steady state calculation (rad  ) [used only when CalcStdy=True]\n');
    fprintf(fid, '   0.0010   VelTol      - Convergence tolerance for the 2-norm of velocities    in the periodic steady state calculation (rad/s) [used only when CalcStdy=True]\n');
    fprintf(fid, '---------------------- MODEL LINEARIZATION -------------------------------------\n');
    fprintf(fid, '  36        NAzimStep   - Number of equally-spaced azimuth steps in periodic linearized model (-)\n');
    fprintf(fid, '   1        MdlOrder    - Order of output linearized model {1: 1st order A, B, Bd, C, D, Dd; 2: 2nd order M, C, K, F, Fd, VelC, DspC, D, Dd} (switch)\n');
    fprintf(fid, '---------------------- INPUTS AND DISTURBANCES ---------------------------------\n');
    fprintf(fid, '   2        NInputs     - Number of control inputs [0 (none) or 1 to 4+NumBl] (-)\n');
    fprintf(fid, '   3,4      CntrlInpt   - List of control inputs [1 to NInputs] {1: nacelle yaw angle, 2: nacelle yaw rate, 3: generator torque, 4: collective blade pitch, 5: individual pitch of blade 1, 6: individual pitch of blade 2, 7: individual pitch of blade 3 [unavailable for 2-bladed turbines]} (-) [unused if NInputs=0]\n');
    fprintf(fid, '   1        NDisturbs   - Number of wind disturbances [0 (none) or 1 to 7] (-)\n');
    fprintf(fid, '   1        Disturbnc   - List of input wind disturbances [1 to NDisturbs] {1: horizontal hub-height wind speed, 2: horizontal wind direction, 3: vertical wind speed, 4: horizontal wind shear, 5: vertical power law wind shear, 6: linear vertical wind shear, 7: horizontal hub-height wind gust} (-) [unused if NDisturbs=0]\n');
    fclose(fid);

    % Wind input file
    fid = fopen('subfunctions\steadywind.hh', 'wt');
    fprintf(fid, '! Created %s.\n', datestr(now));
    fprintf(fid, '! Time Wind  Wind Vert. Horiz. Vert. LinV  Gust\n');
    fprintf(fid, '!      Speed Dir  Speed Shear  Shear Shear Speed\n');
    for t = 0:1
        fprintf(fid, '%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n', ...
            t, WindSpeeds(j), 0, 0, 0, 0, 0, 0);
    end
    fprintf(fid, 'end\n');
    fclose(fid);

    % Run FAST
    system(['"', pwd, '\subfunctions\FAST" [/h] "', pwd, '\subfunctions\linearization.fst"']);

    % Save state-space model
    Name = {'subfunctions\linearization'};
    GetMatsFloat;
    Inputs = [DescCntrlInpt; DescDisturbnc];
    Outputs = [DescOutput];
    States = [DescStates];
    sysm{j} = ss(AvgAMat,[AvgBMat AvgBdMat],AvgCMat,[AvgDMat AvgDdMat],'InputName',Inputs,'Outputname',Outputs);
    LinTorque = str2num(DescCntrlInpt{1}(53:66));
    LinPitch = str2num(DescCntrlInpt{2}(53:66)); 
    Lin.Torque(j) = LinTorque;
    Lin.Pitch(j) = LinPitch;
    Lin.RSpeed(j) = RotSpeed;
    Lin.V(j) = WindSpeeds(j);
        
end

save(get(handles.LinearFile_textbox, 'String'), 'sysm','Lin');
disp(' ')
disp(' ')
disp(' ... Linearization complete!')

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
