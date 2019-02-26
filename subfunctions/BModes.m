function [y11_shape, y11_coeff, y11_freq, ...
    y12_shape, y12_coeff, y12_freq, ...
    y21_shape, y21_coeff, y21_freq, ...
    y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,type,rpm)

%% Blade input files
if type == 1
    
    % Blade input file
    props = [pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'BModes_blade.dat'];
    n_secs = length(Blade.Radius);
    fid = fopen(props, 'wt');
    fprintf(fid, 'Blade section properties\n');
    fprintf(fid, '%0.0f\tn_secs:\tnumber of blade or tower sections at which properties are specified (-)\n\n', n_secs);
    fprintf(fid, 'sec_loc\tstr_tw\ttw_iner\tmass_den\tflp_iner\tedge_iner\tflp_stff\tedge_stff\ttor_stff\taxial_stff\tcg_offst\tsc_offst\ttc_offst\n');
    fprintf(fid, '(-)\t(deg)\t(deg)\t(kg/m)\t\t(kg-m)\t\t(kg-m)\t\t(Nm^2)\t\t(Nm^2)\t\t(Nm^2)\t\t(N)\t\t(m)\t\t(m)\t\t(m)\n');
    for i = 1:n_secs
        fprintf(fid, '%0.5f\t%7.3f\t%7.3f\t%7.2f\t\t%8.2f\t%8.2f\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\n', ...
            (Blade.Radius(i)-Blade.Radius(1))/(Blade.Radius(end)-Blade.Radius(1)), ...
            Blade.Twist(i), ...
            Blade.Twist(i), ...
            Blade.Mass(i), ...
            Blade.FlapIner(i), ...
            Blade.EdgeIner(i), ...
            Blade.EIflap(i), ...
            Blade.EIedge(i), ...
            Blade.GJ(i), ...
            Blade.EA(i), ...
            Blade.cg(i), ...
            Blade.sc(i), ...
            0);
    end
    fclose(fid);
    
end

%% Tower input file
if type == 2
    
    % Tower input file
    props = [pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'BModes_tower.dat'];
    n_secs = length(Tower.Height);
    fid = fopen(props, 'wt');
    fprintf(fid, 'Tower section properties\n');
    fprintf(fid, '%0.0f\tn_secs:\tnumber of blade or tower sections at which properties are specified (-)\n\n', n_secs);
    fprintf(fid, 'sec_loc\tstr_tw\ttw_iner\tmass_den\tflp_iner\tedge_iner\tflp_stff\tedge_stff\ttor_stff\taxial_stff\tcg_offst\tsc_offst\ttc_offst\n');
    fprintf(fid, '(-)\t(deg)\t(deg)\t(kg/m)\t\t(kg-m)\t\t(kg-m)\t\t(Nm^2)\t\t(Nm^2)\t\t(Nm^2)\t\t(N)\t\t(m)\t\t(m)\t\t(m)\n');
    for i = 1:n_secs
        fprintf(fid, '%0.5f\t%0.0f\t%0.0f\t%7.2f\t\t%8.2f\t%8.2f\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\t%1.2E\n', ...
            Tower.Height(i)/Tower.Height(end), ...
            0, ...
            0, ...
            Tower.Mass(i), ...
            Tower.Iner(i), ...
            Tower.Iner(i), ...
            Tower.EI(i), ...
            Tower.EI(i), ...
            Tower.GJ(i), ...
            Tower.EA(i), ...
            0, ...
            0, ...
            0);
    end
    fclose(fid);

end

%% BModes input file
fid = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'BModes.bmi'], 'wt');
fprintf(fid, '======================   BModes v1.03 Main Input File  ==================\n');
fprintf(fid, 'Created %s.\n\n', datestr(now));
fprintf(fid, '--------- General parameters ---------------------------------------------------------------------\n');
fprintf(fid, 'False    Echo        Echo input file contents to *.echo file if true.\n');
fprintf(fid, '%i       beam_type   1: blade, 2: tower (-)\n', type);
fprintf(fid, '%4.1f    romg:       rotor speed, automatically set to zero for tower modal analysis (rpm)\n', rpm);
fprintf(fid, '1.0      romg_mult:  rotor speed muliplicative factor (-)\n');
if type == 1
    fprintf(fid, '%4.2f    radius:     rotor tip radius measured along coned blade axis OR tower height (m)\n', Blade.Radius(end));
    fprintf(fid, '%4.2f    hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\n', Blade.Radius(1));
elseif type == 2
    fprintf(fid, '%4.2f    radius:     rotor tip radius measured along coned blade axis OR tower height (m)\n', Tower.Height(end));
    fprintf(fid, '%4.2f    hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\n', Tower.Height(1));
end
fprintf(fid, '%4.2f    precone:    built-in precone angle, automatically set to zero for a tower (deg)\n', Blade.Cone);
fprintf(fid, '%4.3f    bl_thp:     blade pitch setting, automatically set to zero for a tower (deg)\n', Control.Pitch.Fine);
fprintf(fid, '1         hub_conn:   hub-to-blade connection [1: cantilevered; other options not yet available] (-)\n');
fprintf(fid, '20        modepr:     number of modes to be printed (-)\n');
fprintf(fid, 't         TabDelim    (true: tab-delimited output tables; false: space-delimited tables)\n');
fprintf(fid, 't         mid_node_tw  (true: output twist at mid-node of elements; false: no mid-node outputs)\n\n');
fprintf(fid, '--------- Blade-tip or tower-top mass properties --------------------------------------------\n');
if type == 1
    fprintf(fid, '%4.2f    tip_mass    blade-tip or tower-top mass (kg)\n', 0);
elseif type == 2
    fprintf(fid, '%4.2f    tip_mass    blade-tip or tower-top mass (kg)\n', Nacelle.Housing.Mass + Nacelle.Hub.Mass + Blade.Number*trapz(Blade.Radius,Blade.Mass));
end
fprintf(fid, '0.        cm_loc      tip-mass c.m. offset from the blade axis measured along the tip section y reference axis (m)\n');
fprintf(fid, '0.        ixx_tip     blade lag mass moment of inertia about the tip-section x reference axis (kg-m^2)\n');
fprintf(fid, '0.        iyy_tip     blade flap mass moment of inertia about the tip-section y reference axis (kg-m^2)\n');
fprintf(fid, '0.        izz_tip     torsion mass moment of inertia about the tip-section z reference axis (kg-m^2)\n');
fprintf(fid, '0.        ixy_tip     cross product of inertia about x and y reference axes(kg-m^2)\n');
fprintf(fid, '0.        izx_tip     cross product of inertia about z and x reference axes(kg-m^2)\n');
fprintf(fid, '0.        iyz_tip     cross product of inertia about y and z reference axes(kg-m^2)\n\n');
fprintf(fid, '--------- Distributed-property identifiers --------------------------------------------------------\n');
fprintf(fid, '1         id_mat:     material_type [1: isotropic; non-isotropic composites option not yet available]\n');
fprintf(fid, '"%s" sec_props_file   name of beam section properties file (-)\n\n', props);
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

%% Run BModes
system(['subfunctions' filesep 'bmodes "', pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'BModes.bmi"']);

%% Find modes
% Read output
fileID = fopen([pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'BModes.out'],'r');
textscan(fileID, '%[^\n\r]', 5, 'ReturnOnError', false);
data = textscan(fileID, '%s%s%s%s%s%s%[^\n\r]', 'Delimiter', '\t', 'ReturnOnError', false);
fclose(fileID);
n = (length(data{1}) - 1)/20;
mode = zeros(20,1);
freq = zeros(20,1);
y1 = zeros(20,n-2); % Blade flapwise/tower side-side modes
y2 = zeros(20,n-2); % Blade edgewise/tower fore-aft modes
y3 = zeros(20,n-2); % Torsional modes
for i = 1:20
    header = cell2mat(data{1}((i-1)*n+1));
    freq(i) = str2double(header(end-15:end-4));
    for j = 1:n-2
        y1(i,j) = str2double(data{2}((i-1)*n+2+j));
        y2(i,j) = str2double(data{4}((i-1)*n+2+j));
        y3(i,j) = pi/180*str2double(data{6}((i-1)*n+2+j));
    end
    if max(abs(y1(i,:))) > max(abs(y2(i,:))) && max(abs(y1(i,:))) > max(abs(y3(i,:)))
        mode(i) = 1;
    elseif max(abs(y2(i,:))) > max(abs(y1(i,:))) && max(abs(y2(i,:))) > max(abs(y3(i,:)))
        mode(i) = 2;
    else
        mode(i) = 3;
    end
end
y1 = y1(:,1:2:end);
y2 = y2(:,1:2:end);
y3 = y3(:,1:2:end);
    
% Blade flapwise/tower side-side modes
i = find(mode == 1);
x = linspace(0,1,(n-3)/2+1);
y11_shape = [y1(i(1),:); y2(i(1),:); y3(i(1),:)];
y11_coeff = polyfit(x,y1(i(1),:),6);
y11_coeff = y11_coeff / sum(y11_coeff(1:5));
y11_freq = freq(i(1));
y12_shape = [y1(i(2),:); y2(i(2),:); y3(i(2),:)];
y12_coeff = polyfit(x,y1(i(2),:),6);
y12_coeff = y12_coeff / sum(y12_coeff(1:5));
y12_freq = freq(i(2));

% Blade edgewise/tower fore-aft modes
i = find(mode == 2);
x = linspace(0,1,(n-3)/2+1);
y21_shape = [y1(i(1),:); y2(i(1),:); y3(i(1),:)];
y21_coeff = polyfit(x,y2(i(1),:),6);
y21_coeff = y21_coeff / sum(y21_coeff(1:5));
y21_freq = freq(i(1));
y22_shape = [y1(i(2),:); y2(i(2),:); y3(i(2),:)];
y22_coeff = polyfit(x,y2(i(2),:),6);
y22_coeff = y22_coeff / sum(y22_coeff(1:5));
y22_freq = freq(i(2));
