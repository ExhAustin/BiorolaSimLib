% Plot a drawing of system configuration

function drawFcn(t, x)
% Constants
%***define constants here if needed***
I1 = 0;
I2 = 0;
m1 = 1;
m2 = 1;
L1 = 0.5;
L2 = 0.5;
g = 9.8;

% Plot
%***plot system drawing using x (and t) here***
theta1 = x(1);
theta2 = x(2);

x0 = 0;
y0 = 0;
x1 = x0 + L1*sin(theta1);
y1 = y0 - L1*cos(theta1);
x2 = x1 + L2*sin(theta1+theta2);
y2 = y1 - L2*cos(theta1+theta2);


plot([x0, x1], [y0, y1], 'linewidth', 3);
hold on;
plot([x1, x2], [y1, y2], 'linewidth', 3);

% Resize to fit
center = [0, -0.5]; %***coordinate of the center of the plot (1-by-2 vector)***
width = 2; %***width of the plot window***

screenratio = 768/1366;
x_min = center(1) - width/2;
x_max = center(1) + width/2;
y_min = center(2) - screenratio*width/2;
y_max = center(2) + screenratio*width/2;
axis([x_min, x_max, y_min, y_max]);

end