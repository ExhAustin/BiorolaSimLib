% ode_simulation.m in a function form (returns result structure instead of saving to file)

function result = ode_simulation_func(sim_time, x0);
%%=====================================
% Initialization
%======================================
% Setup state equations & event function
addpath('functions');

odeFcn = @pendulum_sys; %***initial state equation function handle***

eventFcn = @eventFcn_default; %***initial event function handle (@eventFcn_default if no events needed)***
options = odeset('Events', eventFcn);	% add additional ode45 options here if needed

% Time
t_start = 0;
t_finish = sim_time;

% Result containers
T_result = [];
X_result = [];
TE_result = [];
XE_result = [];
IE_result = [];

%%=====================================
% Simulation
%======================================
% Iteratively solve ode
x_start = x0;

while 1
	% Run one instance of ode45
	[T X TE XE IE] = ode45(odeFcn, [t_start t_finish], x_start, options);

	% Record state data
	T_result = [T_result; T];
	X_result = [X_result; X];

	% Record event data (note: events that happen at t_start will be recorded but will not cause ode45 to terminate)
	if length(TE_result) == 0
		idx = 1;
	else
		for i = 1:1:length(TE)
			if TE ~= t_start
				idx = i;
				break;
			end
		end
	end
	TE_result = [TE_result; TE(idx:end)];
	XE_result = [XE_result; XE(idx:end, :)];
	IE_result = [IE_result; IE(idx:end)];


	% Phase transition
		%***change odeFcn (and eventFcn) here according to IE*** (leave blank if no transitions are required)

	% Other
		%***add additional custom settings here*** (ex: reset states, terminate simulation early)

	% End loop if sim_time is reached, otherwise update starting time
	if T(end) ~= sim_time
		% Set new ode initial conditions
		t_start = TE(end);
		x_start = XE(end, :);
	else
		% Simulation complete!
		break;
	end

end

%%=====================================
% Process results
%======================================
% Pack data into a structure
result.dof = length(x0)/2;
result.x0 = x0;
result.T = T_result;
result.X = X_result;
result.TE = TE_result;
result.XE = XE_result;
result.IE = IE_result;
result.datetime = datestr(now, 'yyyy/mm/dd HH:MM:SS');
	%***Add more fields of information about this simulation here(optional)***
end