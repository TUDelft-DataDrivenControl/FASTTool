# FASTTool
The FASTTool has been developed for educational purposes in wind turbine design. FASTTool is a graphical user interface (GUI) for NREL’s aeroelastic simulation code FAST. The tool is centered around a three-dimensional animated wind turbine plot, which dynamically adapts to the defined design inputs. FASTTool provides users with convenient and insightful tools to tune controllers and assess the performance of the design.

The simulation back-end is based on NREL’s FAST v8.16, and is to date still under active development, and updated regularly based on new insights and feedback from students. FASTTool is developed in MATLAB/Simulink, in conjunction with the publicly available FAST simulation code. The choice of software is convenient for an academic environment, since no license fees are demanded for the use of FAST, while often students have access to and experience with MATLAB/Simulink. Although MATLAB and Simulink require a license, the employed environment provides flexibility and insight for both the end-user and developer. The MATLAB scripts and the Simulink model of the tool can be edited by more expert users to add or change functionality and to enable other inputs and outputs.

## Key features
FASTTool provides the following key features:
* Structural, drivetrain and controller design: Blade design allows to to radially specify the blade geometryand structural properties by defining the chord, twist and airfoil for each node, as well as the mass density, flap- and edgewise stiffness. The user can also edit airfoil properties or add new airfoils. A similar interface is provided for tower design. Parameters size the nacelle, and define the drivetrain by eﬃciencies, the gearbox ratio, and the generator inertia. The controller design interface allows to visually tune the pitch controller by loop shaping the system’s frequency responses. 
* Steady-state rotor performance and modal analysis: Modal analyses can be performed on the tower fore-aft and side-side modes, along with the blade flap- and edgewise modes, and visualized using a Campbell diagram.
* Linearization: The non-linear FAST model can be linearized at operating points of interest, and provides functionality for finding trim conditions.
* Wind load cases and simulation: Run certification simulations using various wind conditions, such as steady wind, stepped wind speed changes, a normal or extreme turbulence model.
* Simulink-based controller and simulation environment: During a simulation run, calling FAST dynamic library by a Simulink S-Function, a controller that is provided by Simulink blocks is used, and is configured with information from the different interfaces.

## Referencing
When you use FASTTool in any publication, please cite the following paper:
* Mulders, S.P. and Zaaijer, M.B. and Bos, R. and van Wingerden, J.W. "Wind turbine control: open-source software for control education, standardization and compilation". Journal of Physics: Conference Series. Vol. 1452. No. 1. IOP Publishing, 2020. [Link to the paper](https://iopscience.iop.org/article/10.1088/1742-6596/1452/1/012010)
