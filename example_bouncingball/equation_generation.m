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
m = 1;
g = 9.8;
r = 0.05;

% Generalized coordinates
dof = 2; %***number of degrees of freedom of the system***
syms t real;
q = sym('q', [dof,1]);
q_dot = sym('q_dot', [dof,1]);
x = [q; q_dot];
	
% Names of generalized coordinates
	%***assign names to symbolic variables to increase the readability of the Lagrangian here(optional)***

%%=====================================
% Generate state equation files
%======================================
	% Note: copy & paste several versions of the following code to generate multiple state equations
%-------------------------------------------------------------
% Lagrangian
T = 0.5*m*(q_dot(1)^2 + q_dot(2)^2);
V = m*g*q(2);
L = T - V; %***Lagrangian of the system here***

% Generalized force
tau = 0; %***generalized force here*** (as a column vector of functions of t, q, dq)
	% Note: tau = 0 for passive systems

% Constraints
constraints = 0; %***constraints here*** (as a column vector of functions of t, q)
	% Note: constraints = 0 for unconstrained systems

% Calculate Lagrangian (also checks for errors)
f = lagrangian2eqns(t, q, q_dot, L, tau, constraints);

% Save symbolic function as a Matlab function file
filename = 'ball_sys'; %***name of generated function here***
pathname = strcat('functions/', filename);
matlabFunction(f, 'file', pathname, 'vars', {t, x});
%-------------------------------------------------------------

%%=====================================
% Generate event function files
%======================================
	% Note1: copy & paste several versions of the following code to generate multiple event functions
	% Note2: Remove a "%" from the next line to comment this entire section if no event functions are needed
%%{
%-------------------------------------------------------------
% Setup
n_events = 1; %***number of events the system might encounter at this phase***

value = q(2) - r; %***value of event function (n_events-by-1 column vector of functions of q and/or q_dot)***
isterminal = 1; %***if the integration should terminate at a zero of this event function (1:yes, 2:no) (n_events-by-1 column vector)***
direction = -1; %***direction of value at a zero (1: value is increasing, -1: value is decreasing, 0: both) (n_events-by-1 column vector)***

% Save symbolic function as a Matlab function file
filename = 'ground_event'; %***name of generated function here***
pathname = strcat('functions/', filename);
matlabFunction(value, isterminal, direction, 'file', pathname, 'vars', {t, x});
%-------------------------------------------------------------
%}