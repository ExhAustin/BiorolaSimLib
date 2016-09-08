% Plot a drawing of system configuration

function drawFcn(t, x)
% Constants
%***define constants here if needed***

% Plot
%***plot system drawing using x (and t) here***

% Resize to fit
center = %***coordinate of the center of the plot (1-by-2 vector)***
width = %***width of the plot window***

screenratio = 768/1366;
x_min = center(1) - width/2;
x_max = center(1) + width/2;
y_min = center(2) - screenratio*width/2;
y_max = center(2) + screenratio*width/2;
axis([x_min, x_max, y_min, y_max]);

end