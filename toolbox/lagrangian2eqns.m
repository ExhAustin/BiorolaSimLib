%=========================================================================
% lagrangian2eqns.m
%=========================================================================
% Author: Austin Wang
% Year: 2016
% Requirements: MATLAB version R2013a or newer, Symbolic Math Toolbox
%
% Description:
% 	Derives state equations of a system from its Lagrangian.
%	Requires MATLAB Symbolic Math Toolbox.
%
% Syntax:
%	f = lagrangian2eqns(t, q, q_dot, L)
%	f = lagrangian2eqns(t, q, q_dot, L, tau)
%	f = lagrangian2eqns(t, q, q_dot, L, tau, constr)
%
% Arguments:
%	t
%		Time variable.
%		Specified as a single symbolic variable
%
%	q
%		Generalized coordinates of the system.
%		Specified as a n_dof-by-1 column vector of symbolic variables.	(n_dof = number of dof of the system)
%
%	q_dot
%		Generalized velocities of the system.
%		Specified as a n_dof-by-1 column vector of symbolic variables.
%
%	L
%		Lagrangian of the system.
%		Specified as a symbolic expression constructed by symbolic variables in q and q_dot.
%
%	tau
%		Generalized force.
%		If no force is applied, specified as 0,
%		else specified as a n_dof-by-1 column vector of symbolic expressions constructed by t and/or symbolic variables in q, and q_dot.
%
%	constr
%		Holonomic constraints in the form of constr==0.
%		If no constraints exist, specified as 0,	
%		else specified as a n-by-1 column vector of symbolic expressions constructed by t and/or symbolic variables in q.
%
%	f
%		The system state equations as a symbolic function. 
%			dx/dt = f(x), where x = [q; q_dot]
%
%
% Example: Point mass in gravity
%
%	m = 1;
%	g = 9.8;
%
%	syms t real;
%	
%	dof = 2;
%	q = sym('q', [dof,1]);		% q(1) = x, q(2) = y
%	q_dot = sym('q_dot', [dof,1]);
%	
%	T = 0.5*m*( q_dot(1)^2 + q_dot(2)^2 );
%	V = m*g*q(2);
%	L = T - V;
%	
%	f = lagrangian2eqns(t, q, q_dot, L);
%	
%	x = [q; q_dot];
%	matlabFunction(f, 'file', '2dpointmass_sys', 'vars', {t, x});
%
%=========================================================================
function f = lagrangian2eqns(t, q, q_dot, L, varargin)
% Parse input
[n_dof, tau, constraints, n_constraints] = parseInput(t, q, q_dot, L, varargin, nargin);

% Formulate manipulator equations ( H(q)*ddq + C(dq,q)*dq + G(q) == tau )
H = jacobian(jacobian(L, q_dot), q_dot);
C = jacobian(jacobian(L, q_dot), q);
G = - jacobian(L, q).';

% Solve Lagrange multipliers of holonomic constraints
if n_constraints > 0
	lambda = sym('lambda', [n_constraints, 1]);
	q_ddot = sym('q_ddot', [n_dof, 1]);
	vars = num2cell([lambda; q_ddot].');

	sys_eqns = simplify( H*q_ddot + C*q_dot + G - jacobian(constraints, q).'*lambda - tau );

	q_str = q;
	constraints_str = constraints;
	for i = 1:n_dof
		q_str(i) = strcat('q',num2str(i),'(t)');
		constraints_str = subs(constraints_str, q(i), q_str(i));
	end
	constraint_eqns = diff(diff(constraints_str, t), t);
	for i = 1:n_dof
		constraint_eqns = subs(constraint_eqns, diff(q_str(i),t,2), q_ddot(i));
		constraint_eqns = subs(constraint_eqns, diff(q_str(i),t), q_dot(i));
		constraint_eqns = subs(constraint_eqns, q_str(i), q(i));
	end
	constraint_eqns = simplify(constraint_eqns);

	S = solve([sys_eqns; constraint_eqns], vars{:});
	S_cell = struct2cell(S);
	for i = 1:n_constraints
		lambda_ans(i,1) = S_cell{i};
	end

	tau_constraint = simplify( jacobian(constraints, q).'*lambda_ans );
else
	tau_constraint = zeros(n_dof, 1);
end

% Calculate system equations using jacobians ( dx/dt = f(t, x) )
f = simplify([
	q_dot;
	H \ ( tau + tau_constraint - G - C*q_dot );
]);

end

% Function for checking input validity
function [n_dof, tau, constraints, n_constraints] = parseInput(t, q, q_dot, L, vararg, n_arg)	
	% Parse variable length input arguments
	if n_arg > 6
		error('Error: Too many input arguments.');
	elseif n_arg < 4
		error('Error: Not enough input arguments.')
	elseif n_arg == 4
		tau = 0;
		constraints = 0;
	elseif n_arg == 5
		tau = vararg{1};
		constraints = 0;
	elseif n_arg == 6 
		tau = vararg{1};
		constraints = vararg{2};
	end

	% t
	if ~isequal(size(t),[1,1]) || ~strcmp(class(t), 'sym')
		error('Error: The time variable t must be a single symbolic variable.');
	end

	% q
	n_dof = length(q);
	if ~isequal(size(q),[n_dof,1]) || ~strcmp(class(q), 'sym')
		error('Error: The generalized coordinates q must be a n_dof-by-1 column vector of symbolic variables.');
	end

	% q_dot
	if ~isequal(size(q_dot),[n_dof,1]) || ~strcmp(class(q_dot), 'sym')
		error('Error: The derivative of generalized coordinates q_dot must be a n_dof-by-1 column vector of symbolic variables.');
	end

	% L
	if ~isequal(size(L),[1,1]) || ~strcmp(class(L), 'sym')
		error('Error: The Lagrangian L must be a single symbolic expression.');
	end
	
	symvec_L = symvar(L);
	if ~all( ismember(symvec_L, [q; q_dot]) )
		error('Error: The Lagrangian L must be constructed by symbolic variables in q and q_dot.')
	end
	
	% tau
	if ~isequal(size(tau), [n_dof, 1])
		if ~( isequal(size(tau), [1, 1]) && tau == 0 )
			error('Error: The generalized force tau must be 0 or a n_dof-by-1 column vector.');
		else
			tau = zeros(n_dof, 1);
		end
	else
		symvec_tau = symvar(tau);
		if ~all( ismember(symvec_tau, [t; q; q_dot]) )
			error('Error: The generalized force tau must be constructed by t and/or symbolic variables in q and q_dot.')
		end
	end

	% constraints
	symvec_constraints = symvar(constraints);
	if isequal(size(constraints), [1, 1]) && constraints == 0
		n_constraints = 0;
	elseif size(constraints,2) ~= 1 || ~strcmp(class(L), 'sym')
		error('Error: Constraints must be 0 or a n-by-1 column vector of symbolic expressions.')
	elseif ~all( ismember(symvec_constraints, [t; q]) )
		error('Error: Constraints must be constructed by t and/or symbolic variables in q.')
	else
		n_constraints = length(constraints);
	end
end