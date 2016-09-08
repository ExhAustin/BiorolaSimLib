%=========================================================================
% downhillOpt.m
%=========================================================================
% Author: Austin Wang
% Year: 2016
% Requirements: MATLAB version R2013a or newer
%
% Description:
% 	Optimizes a function by recursively searching for better results inside a user-defined mesh.
%	The algorithm starts an initial point, then repeatly picks the neighboring node with lowest cost and continues in that direction.
%	Multiple initial points can be chosen.
%	If the cost is lower than the threshold, the algorithm then uses DFS to find all nodes with costs lower than the threshold in the area.
%
% Syntax:
% 	[x_opt, cost_opt, y_opt] = downhillOpt(fun, x0, Mesh)
% 	[x_opt, cost_opt, y_opt] = downhillOpt(fun, x0, Mesh, threshold)
%
% Arguments:
%	fun
%		The function to be minimized.
% 		It accepts input x and returns a cost as its first output argument.
%		If the function fun has more than 1 output arguments, the 2nd output will be recorded in y_opt.
%		Specified as a function handle for a function file (ex: @myfun ),
%		or as a function handle for an anonymous function (ex: @(x,y)sin(x^2+y^2) ).
%		Note: This algorithm will not be effective if fun is a chaotic function.
%
%	x0
%		The initial starting point of the search.
%		Can be specified as a vector (single initial point) or a cell array of vectors (multiple initial points).
%
%	Mesh
%		All possibilites of function inputs, which can be thought of a mesh in n-dimensional space,
%		with each dimension as an input parameter.
%		Specified as a cell array of vectors.
%
%	threshold
%		If specified as an real number, downhillOpt will return all function inputs that result in a cost lower than the threshold.
%		If not specified, downhillOpt will find the function input that results in the lowest cost.
%		Note: If no cost lower than the threshold is found, downhillOpt returns the input that results in the lowest cost it can find.
%		Note2: A higher threshold will have a better chance to find the optimal result(s), while increasing program runtime.
%
%	x_opt
%		The optimal function input(s) found.
%		A cell array of vectors.
%
%	cost_opt
%		The optimal cost found.
%		A cell array of real numbers.
%
%	y_opt
%		The secondary function output(s) produced by the optimal function input(s) found.
%		A cell array of the output data type of fun.
%
% Example:
%
%	fun = @(x,y)sin(x^2+y^2);
%	x0 = {
%		[0, 0],
%		[1.5, -1.5]
%	};
%	Mesh = {
%		[-4:0.001:4],
%		[-4:0.001:4]
%	};
%	threshold = -0.5;
% 	[x_opt, cost_opt, y_opt] = downhillOpt(fun, x0, Mesh, threshold);
%
%=========================================================================
function [x_opt, varargout] = downhillOpt(fun, x0, Mesh, varargin)
	%%=====================================
	% Initialization
	%======================================
	% Parse input
	[n_dim, fun_nargout, findmin, threshold, multiple_x0] = parseInput(fun, x0, Mesh, varargin, nargin);
	
	% Mesh dimensions
	for i = 1:n_dim
		mesh_size(i) = numel(Mesh{i});
	end
	
	mesh_size_list = num2cell(mesh_size);
	
	% Initialize cost database
	CostData = zeros(mesh_size_list{:});
	
	% Initialize visit record database
	VisitData = zeros(mesh_size_list{:});
			% 0: unevaluated, 1: evaluated, 2: visited & evaluated
	
	% Initialize stack
	SearchStack = [];
	
	if multiple_x0
		for i = 1:numel(x0)
			if size(x0{i},1) <= size(x0{i},2)
				SearchStack = [SearchStack ; x2coord( x0{i} ) ];
			else
				SearchStack = [SearchStack ; x2coord( transpose(x0{i}) ) ];
			end
		end
	else
		if size(x0,1) <= size(x0,2)
			SearchStack = x2coord( x0 );
		else
			SearchStack = x2coord( transpose(x0) );
		end
	end
	
	% Initialize counters
	n_pass = 0;
	n_steps = 0;
	n_evaluations = 0;
	
	% Initialize minimum cost
	min_cost = inf;
	
	% Execution time measurement
	start_time = clock;
	
	%%=====================================
	% Downhill Iteration
	%======================================
	clc;
	fprintf('Iterating to local minima...\n');
	
	while isempty(SearchStack) ~= 1
		current = SearchStack(end,:);
		SearchStack(end,:) = [];				% SearchStack.pop
	
		% Check if this coordinate has been visited before
		current_list = num2cell(current);
		if VisitData(current_list{:}) == 2
			continue;
		end
	
		% Deviate each parameter coordinate
		TestQueue = current;
		for i = 1:1:n_dim			% alters ith parameter
			TestQueue_new = [];
			while ~isempty(TestQueue)
				target_old = TestQueue(end,:);
				TestQueue(end,:) = [];									%TestQueue.dequeue
				for deviation = [1 -1 0]		
					if target_old(i)+deviation > 0 && target_old(i)+deviation < mesh_size(i)+1		%need if new coordinate is within range
						target_new = target_old;
						target_new(i) = target_old(i) + deviation;
						TestQueue_new = [target_new ; TestQueue_new];	%TestQueue_new.enqueue
					end
				end
			end
			TestQueue = TestQueue_new;
		end
	
		% Obtain function evaluation results of all coordinates in TestQueue
		n_neighbors = size(TestQueue, 1);
		for i_test = 1:n_neighbors
			coordinate_list = num2cell(TestQueue(i_test,:));
			if VisitData(coordinate_list{:}) == 0
				% Evaluate target function
				for i = 1:1:n_dim
					x(1,i) = Mesh{i}(TestQueue(i_test,i));
				end
				x_list = num2cell(x);

				if fun_nargout > 1
					[cost, y] = fun(x_list{:});
				else
					cost = fun(x_list{:});
				end
				n_evaluations = n_evaluations + 1;

				group_cost(i_test) = cost;
				CostData(coordinate_list{:}) = cost;
	
				% Record if it is the minimum cost so far
				if cost < min_cost
					min_cost = cost;
	
					if findmin
						x_opt = x;
						cost_opt = cost;
						if fun_nargout > 1
							y_opt = y;
						end
					elseif n_pass == 0;
						x_opt{1} = x;
						cost_opt{1} = cost;
						if fun_nargout > 1
							y_opt{1} = y;
						end
					end
				end
	
				% Record if cost is lower than threshold
				if ~findmin
					if cost <= threshold
						n_pass = n_pass + 1;
	
						x_opt{n_pass} = x;
						cost_opt{n_pass} = cost;
						if fun_nargout > 1
							y_opt{n_pass} = y;
						end
					end
				end
	
				% Mark as evaluated
				VisitData(coordinate_list{:}) = 1;
			else
				% Obtain cost from database
				group_cost(i_test) = CostData(coordinate_list{:});
			end
		end
	
		% Mark current as visited
		VisitData(current_list{:}) = 2;
		n_steps = n_steps + 1;

		% for every neighbour:
		% if pass && not visited => push into SearchStack
		all_fail = 1;
		for i_test = 1:1:n_neighbors
			if group_cost(i_test) <= threshold
				all_fail = 0;
	
				neighborCoord_list = num2cell(TestQueue(i_test,:));
				if VisitData(neighborCoord_list{:}) ~= 2
					SearchStack = [SearchStack ; TestQueue(i_test,:)];			%SearchStack.push
				end
			end
		end
		% if all fail => push neighbour with lowest cost into SearchStack (if not visited)
		if all_fail == 1
			[dummy max_direction] = min(group_cost);
			direction_list = num2cell(TestQueue(max_direction,:));
			if VisitData(direction_list{:}) ~= 2
				SearchStack = [SearchStack ; TestQueue(max_direction,:)];		%SearchStack.push
			end
		end
	
		% Display progress
		clc;
		fprintf('Iterating to local minima...\n');
		fprintf('Function evaluations: %d\n', n_evaluations);
		fprintf('Steps taken: %d\n', n_steps);
		fprintf('Minimum cost: %f\n', min_cost);
		if ~findmin
			fprintf('Successful sets of parameters found: %d\n\n', n_pass);
		end
	
	end
	
	%%=====================================
	% Finalize & Output
	%======================================
	% Report execution time
	finish_time = clock;
	exe_time = etime(finish_time, start_time);
	exe_hr = floor(exe_time/3600);
	exe_time = mod(exe_time, 3600);
	exe_min = floor(exe_time/60);
	exe_time = mod(exe_time, 60);
	fprintf('Execution time: %d hours, %d minutes, %2.2f seconds\n', exe_hr, exe_min, exe_time);
	
	% evaluation result output
	if fun_nargout == 1
		y_opt = {};
	end

	varargout{1} = cost_opt;
	varargout{2} = y_opt;

	%%=====================================
	% Nested functions
	%======================================
	% Locate closest node in mesh
	function coord = x2coord(x_in)
		for i = 1:n_dim
			[dummy, closest_node] = min(abs(Mesh{i} - x_in(i)));
			coord(i) = closest_node;
		end
	end
end

%%=====================================
% Input Parser
%======================================
function [n_dim, fun_nargout, findmin, threshold, multiple_x0] = parseInput(fun, x0, Mesh, vararg, n_arg)
	% Parse variable length input arguments
	if n_arg > 4
		error('Error: Too many input arguments.');
	elseif n_arg < 3
		error('Error: Not enough input arguments.')
	elseif n_arg == 3
		findmin = true;
		threshold = -inf;
	elseif n_arg == 4
		findmin = false;
		threshold = vararg{1};
	end
	
	% fun
	if ~strcmp(class(fun), 'function_handle')
		error('Error: The optimization target must be a function handle.');
	end
	
	n_dim = nargin(fun);
	fun_nargout = nargout(fun);
	
	if n_dim == 0
		error('Error: The optimization target function must have input parameters.');
	elseif n_dim < 0
		error('Error: The inputs of the desired optimization target must not have variable-length input argument list "varargin".');
	end
	
	% x0
	if iscell(x0)
		multiple_x0 = true;

		for i = 1:numel(x0)
			if ~isnumeric(x0{i}) || ( size(x0{i},1)~=1 && size(x0{i},2)~=1 )
				error('Error: The initial input set(s) x0 must be a vector or cell array of vectors');
			end
		end
	elseif isnumeric(x0)
		multiple_x0 = false;

		if size(x0,1) ~= 1 && size(x0,2) ~= 1
			error('Error: The initial input set(s) x0 must be a vector or cell array of vectors');
		end
	else
		error('Error: The initial input set(s) x0 must be a vector or cell array of vectors');
	end

	% Mesh
	if ~iscell(Mesh)
		error('Error: The input parameter mesh must be a cell array of vectors.');
	else
		if numel(Mesh) ~= n_dim
			error('Error: The mesh must have the same dimensions as the inputs of the optimization target function.');
		end

		for i = 1:numel(Mesh)
			if ~isnumeric(Mesh{i})
				error('Error: The input parameter mesh must be a cell array of vectors.');
			end
		end
	end

	% threshold
	if numel(threshold) ~= 1
		error('Error: The cost function threshold must be a single real number.');
	end

end