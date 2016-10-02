%=========================================================================
% equation_generation.m
%
% Utilizes MATLAB symbolic toolbox to generate state equations of a system
% in the form of a MATLAB function file using the Lagrangian method.
%
% Dependencies:
%	toolbox/lagrangian2eqns.m
%=========================================================================
close all;
clear all;
clc;
%%=====================================
% Initialization
%======================================
% Add toolbox path (lagrangian2eqns)
addpath('toolbox');

%%=====================================
% System characteristics
%======================================
% Constants
%***define constants here if needed***
I1 = 0;
I2 = 0;
m1 = 1;
m2 = 1;
L1 = 0.5;
L2 = 0.5;
g = 9.8;

% Generalized coordinates
dof = 2; %***number of degrees of freedom of the system***
syms t real;
q = sym('q', [dof,1]);
q_dot = sym('q_dot', [dof,1]);
x = [q; q_dot];
	
% Names of generalized coordinates
	%***assign names to symbolic variables to increase the readability of the Lagrangian here(optional)***
theta1 = q(1);
theta2 = q(2);
theta1_dot = q_dot(1);
theta2_dot = q_dot(2);

%%=====================================
% Generate state equation files
%======================================
	% Note: copy & paste several versions of the following code to generate multiple state equations
%-------------------------------------------------------------
% Lagrangian
T_trans = 0.5*m1*((L1/2)*theta1_dot)^2 + 0.5*m2*( ...
	( L1*theta1_dot*sin(theta1) + (L2/2)*theta2_dot*sin(theta1+theta2) )^2 ...
	+ ( L1*theta1_dot*cos(theta1) + (L2/2)*theta2_dot*cos(theta1+theta2) )^2 ...
);
T_rot = 0.5*I1*(theta1_dot)^2 + 0.5*I2*(theta1_dot+theta2_dot)^2;
V = -m1*g*(L1/2)*cos(theta1) - m2*g*( (L1)*cos(theta1) + (L2/2)*cos(theta1+theta2) );
L = T_trans + T_rot - V; %***Lagrangian of the system here***

% Generalized force
tau = 0; %***generalized force here*** (as a column vector of functions of t, q, dq)
	% Note: tau = 0 for passive systems

% Constraints
constraints = 0; %***constraints here*** (as a column vector of functions of t, q)
	% Note: constraints = 0 for unconstrained systems

% Calculate Lagrangian (also checks for errors)
f = lagrangian2eqns(t, q, q_dot, L, tau, constraints);

% Save symbolic function as a Matlab function file
filename = 'pendulum_sys'; %***name of generated function here***
pathname = strcat('functions/', filename);
matlabFunction(f, 'file', pathname, 'vars', {t, x});
%-------------------------------------------------------------

%%=====================================
% Generate event function files
%======================================
	% Note1: copy & paste several versions of the following code to generate multiple event functions
	% Note2: Remove a "%" from the next line to comment this entire section if no event functions are needed
%{
%-------------------------------------------------------------
% Setup
n_events = %***number of events the system might encounter at this phase***

value = %***value of event function (n_events-by-1 column vector of functions of q and/or q_dot)***
isterminal = %***if the integration should terminate at a zero of this event function (1:yes, 0:no) (n_events-by-1 column vector)***
direction = %***direction of value at a zero (1: value is increasing, -1: value is decreasing, 0: both) (n_events-by-1 column vector)***

% Save symbolic function as a Matlab function file
filename = %***name of generated function here***
pathname = strcat('functions/', filename);
matlabFunction(value, isterminal, direction, 'file', pathname, 'vars', {t, x});
%-------------------------------------------------------------
%}