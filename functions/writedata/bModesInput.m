function text = bModesInput(Blade, Tower, Nacelle, Control, type, speed, fileLocation)
% bModesInput writes data in a character array, which can be written to
% file. The writing of the character array is separated by this function
% from writing the file in FileIOClass.m, to separate agnostic fileIO from
% the generation of the content of the file.
%
%   Syntax
%     text = bModesInput(Blade, Tower, Nacelle, Control, type, speed, fileLocation)
%
%   Input arguments
%     Blade - Structure with blade specifications
%     Tower - Structure with tower specifications
%     Nacelle - Structure with nacelle specifications, used to determine
%     the tower top mass
%     Control - Structure with control specifications
%     type - Beam-type specification: blade or tower
%     speed - Rotor speed, relevant for centrifugal stiffening (not used
%     for tower)
%     fileLocation - Location of input file for bModes with blade data or
%     tower data, generated using bModesBlade.m or bModesTower.m
%
%   Output
%     text - Character array formatted as an BModes v1.03 main input file
%     that can be written to a plain-text file named BModes.bmi.
%
%   Behaviour
%     Relevant data from the WindTurbineClass properties is printed to a
%     character array according to the specified format. Several parameter
%     values are not taken from the input data, but are hard coded. This is
%     particularly the case for most model parameters or settings. Some
%     parameters are given a constant or calculated value.
%
%   Called by
%     WindTurbineClass.performModalAnalysis
%     WindTurbineClass.linearise
%     WindTurbineClass.simulate

text = '';
text = append(text, sprintf('======================   BModes v1.03 Main Input File  ==================\n'));
text = append(text, sprintf('Created %s.\n\n', datetime));
text = append(text, sprintf('--------- General parameters ---------------------------------------------------------------------\n'));
text = append(text, sprintf('False    Echo        Echo input file contents to *.echo file if true.\n'));
text = append(text, sprintf('%i       beam_type   1: blade, 2: tower (-)\n', type));
text = append(text, sprintf('%4.1f    romg:       rotor speed, automatically set to zero for tower modal analysis (rpm)\n', speed));
text = append(text, sprintf('1.0      romg_mult:  rotor speed muliplicative factor (-)\n'));
if type == 1
    text = append(text, sprintf('%4.2f    radius:     rotor tip radius measured along coned blade axis OR tower height (m)\n', Blade.Radius(end)));
    text = append(text, sprintf('%4.2f    hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\n', Blade.Radius(1)));
elseif type == 2
    text = append(text, sprintf('%4.2f    radius:     rotor tip radius measured along coned blade axis OR tower height (m)\n', Tower.Height(end)));
    text = append(text, sprintf('%4.2f    hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)\n', Tower.Height(1)));
end
text = append(text, sprintf('%4.2f    precone:    built-in precone angle, automatically set to zero for a tower (deg)\n', Blade.Cone));
text = append(text, sprintf('%4.3f    bl_thp:     blade pitch setting, automatically set to zero for a tower (deg)\n', Control.Pitch.Fine));
text = append(text, sprintf('1         hub_conn:   hub-to-blade connection [1: cantilevered; other options not yet available] (-)\n'));
text = append(text, sprintf('20        modepr:     number of modes to be printed (-)\n'));
text = append(text, sprintf('t         TabDelim    (true: tab-delimited output tables; false: space-delimited tables)\n'));
text = append(text, sprintf('t         mid_node_tw  (true: output twist at mid-node of elements; false: no mid-node outputs)\n\n'));
text = append(text, sprintf('--------- Blade-tip or tower-top mass properties --------------------------------------------\n'));
if type == 1
    text = append(text, sprintf('%4.2f    tip_mass    blade-tip or tower-top mass (kg)\n', 0));
elseif type == 2
    text = append(text, sprintf('%4.2f    tip_mass    blade-tip or tower-top mass (kg)\n', Nacelle.Housing.Mass + Nacelle.Hub.Mass + Blade.Number*trapz(Blade.Radius,Blade.Mass)));
end
text = append(text, sprintf('0.        cm_loc      tip-mass c.m. offset from the blade axis measured along the tip section y reference axis (m)\n'));
text = append(text, sprintf('0.        ixx_tip     blade lag mass moment of inertia about the tip-section x reference axis (kg-m^2)\n'));
text = append(text, sprintf('0.        iyy_tip     blade flap mass moment of inertia about the tip-section y reference axis (kg-m^2)\n'));
text = append(text, sprintf('0.        izz_tip     torsion mass moment of inertia about the tip-section z reference axis (kg-m^2)\n'));
text = append(text, sprintf('0.        ixy_tip     cross product of inertia about x and y reference axes(kg-m^2)\n'));
text = append(text, sprintf('0.        izx_tip     cross product of inertia about z and x reference axes(kg-m^2)\n'));
text = append(text, sprintf('0.        iyz_tip     cross product of inertia about y and z reference axes(kg-m^2)\n\n'));
text = append(text, sprintf('--------- Distributed-property identifiers --------------------------------------------------------\n'));
text = append(text, sprintf('1         id_mat:     material_type [1: isotropic; non-isotropic composites option not yet available]\n'));
text = append(text, sprintf('"%s" sec_props_file   name of beam section properties file (-)\n\n', fileLocation));
text = append(text, sprintf('Property scaling factors..............................\n'));
text = append(text, sprintf('1.0       sec_mass_mult:   mass density multiplier (-)\n'));
text = append(text, sprintf('1.0       flp_iner_mult:   blade flap or tower f-a inertia multiplier (-)\n'));
text = append(text, sprintf('1.0       lag_iner_mult:   blade lag or tower s-s inertia multiplier (-)\n'));
text = append(text, sprintf('1.0       flp_stff_mult:   blade flap or tower f-a bending stiffness multiplier (-)\n'));
text = append(text, sprintf('1.0       edge_stff_mult:  blade lag or tower s-s bending stiffness multiplier (-)\n'));
text = append(text, sprintf('1.0       tor_stff_mult:   torsion stiffness multiplier (-)\n'));
text = append(text, sprintf('1.0       axial_stff_mult: axial stiffness multiplier (-)\n'));
text = append(text, sprintf('1.0       cg_offst_mult:   cg offset multiplier (-)\n'));
text = append(text, sprintf('1.0       sc_offst_mult:   shear center multiplier (-)\n'));
text = append(text, sprintf('1.0       tc_offst_mult:   tension center multiplier (-)\n\n'));
text = append(text, sprintf('--------- Finite element discretization --------------------------------------------------\n'));
text = append(text, sprintf('20        nselt:     no of blade or tower elements (-)\n'));
text = append(text, sprintf('Distance of element boundary nodes from blade or flexible-tower root (normalized wrt blade or tower length), el_loc()\n'));
text = append(text, sprintf('0.0	0.05	0.1	0.15	0.2	0.25	0.3	0.35	0.4	0.45	0.5	0.55	0.6	0.65	0.7	0.75	0.8	0.85	0.9	0.95	1.0\n\n'));
text = append(text, sprintf('--------- Properties of tension wires suporting the tower --------------------------------\n'));
text = append(text, sprintf('0         n_attachments: no of wire-attachment locations on tower, maxm allowable is 2; 0: no tension-wire support (-)\n'));
text = append(text, sprintf('3 3       n_wires:       no of wires attached at each location (must be 3 or higher) (-)\n'));
text = append(text, sprintf('6 9       node_attach:   node numbers of attacments location (node number must be more than 1 and less than nselt+2) (-)\n'));
text = append(text, sprintf('0.e0 0.e0 wire_stfness:  wire spring constant in each set (N/m)\n'));
text = append(text, sprintf('0. 0.     th_wire:       angle of tension wires wrt the tower axis at each attachment point (deg)\n'));
