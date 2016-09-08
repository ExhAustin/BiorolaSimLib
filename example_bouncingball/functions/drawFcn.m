% Plot a drawing of system configuration

function drawFcn(t, x)
% Constants
%***define constants here if needed***
m = 1;
g = 9.8;
r = 0.05;

% Plot
%***plot system drawing using x (and t) here***
plot([-5, 5], [0, 0], 'linewidth', 3);
hold on;
viscircles([x(1), x(2)], r);

% Resize to fit
center = [0, 0.5]; %***coordinate of the center of the plot (1-by-2 vector)***
width = 2; %***width of the plot window***

screenratio = 768/1366;
x_min = center(1) - width/2;
x_max = center(1) + width/2;
y_min = center(2) - screenratio*width/2;
y_max = center(2) + screenratio*width/2;
axis([x_min, x_max, y_min, y_max]);

end