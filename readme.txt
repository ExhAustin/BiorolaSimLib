=========¡¸,:*¡E\(^ £s ^)/¡E:*¡E¡¸=========
   Biorola Simulation Library/Template
==========================================
Author: Austin Wang
Year: 2016
Requirements: MATLAB version R2013a or newer, Symbolic Math Toolbox
	(Note: R2014a or newer is preferred due to GUI issues)

=== Introduction ===
A template for deriving the state equations of a system and simulating it through time.
In template files, areas enclosed in "***...***" should be filled in by the user.
Two examples are included in this package. (double pendulum and bouncing ball)
Functionalities in the toolbox can also be explicitly used without following the templates.

=== File List ===
Toolbox files (packaged functions that perform complex tasks):

	toolbox\downhillMinSearch.m
		A downhill search algorithm that finds the local minimum of a function.
		(Specific information is included in the comments in the file)

	toolbox\lagrangian2eqns.m
		Formulates the system state equations from the Lagrangian.
		(Specific information is included in the comments in the file)

	toolbox\plot_gui.m
		A interactive GUI for visualizing simulation results.
		(Generated video files will be placed in \videos)

	toolbox\plot_gui_old.m
		A interactive GUI for visualizing simulation results.
		For R2014a and older.
		(GUI is flawed, can only output to video file)
		(Generated video files will be placed in \videos)

Template files (incomplete codes that users can modify):

	equation_generation.m
		Generates system equations and event functions.
		(Generated files will be placed in \functions)

	ode_simulation.m
		Utilizes ode45 to simulate an ode system.
		(Generated files will be placed in \results)

	functions\drawFcn.m
		Plot a drawing of the system configuration from generalized coordinates.

	functions\ode_simulation_func.m
		Function form of ode_simulation.m (for mass searching).
		Returns simulation results as output instead of saving to file.

Other Files

	result_visualizer.m
		Calls respective GUI functions to visualize simulation results.

	functions\eventFcn_default.m
		An empty event function as a default for systems without events.

=== Instruction Manual ===
1. For each new system, make a copy of the entire \template folder
2. Define system equations in equation_generation.m
3. Setup simulation settings and initial conditions in ode_simulation.m
4. Define a plot of system in functions\drawFcn.m
5. Run equationGeneration.m to generate equation files
6. Run ode_simulation.m to obtain simulation results
7. Run result_visualizer.m to view/save simulation results