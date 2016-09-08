% The default event detection function called by ode45 when no events are to be detected.
function [value isterminal direction] = eventFcn_default(t, x)
	value = 1;
	isterminal = 0;
	direction = 1;
end