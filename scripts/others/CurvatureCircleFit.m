%-------------------------------------------------------------------------------------------------------------------------------
%                                                           OBJECTIVE                                                           
%-------------------------------------------------------------------------------------------------------------------------------
% Here is a brief summary of what the code appears to be doing based on ChatGPT interpretation:
% - The code defines some initial conditions such as a range of values for the Mach number (M), the specific heat ratio (gamma_eos),
% the ambient pressure (pinf), and the freestream velocity (ua).
% - It then calculates the density and velocity of a water droplet at the stagnation point (uw) using the given Mach number and
% assuming incompressible flow.
% - The code defines a range of angles (theta_droplet) and uses them to plot the shape of the water droplet (a circle).
% - The code then proceeds to do some ray tracing calculations to determine the path of light rays passing through the water droplet.
% It appears that the code is iterating over different values of n (the refractive index of water) and alpha (the angle of incidence
% of the light ray) to calculate the position of the light ray after passing through the droplet. The resulting positions are plotted
% in a loop for each value of n.
% - Finally, the code fits a circle to the positions of the light rays using the "circfit" function, which returns the radius of the
% circle (r) and its center (xc, yc). The code then plots the value of r normalized by the droplet radius (11e-3) as a function of
% Mach number.
% It is worth noting that the code is currently commented out in some places (e.g., the last line of the code block).
% Therefore, it is possible that some parts of the code are not currently executing or are incomplete.

%-------------------------------------------------------------------------------------------------------------------------------
%                                                              CODE                                                           
%-------------------------------------------------------------------------------------------------------------------------------
clear all; close all; clc;

%% Initial conditions
M = [1:0.01:10];

gamma_eos = 2.35;       
pinf = 1.00e9; 
ua = 345*M;
p1 = 1.01325e5;
rho1 = 1028;
uw = sqrt(gamma_eos*(pinf+p1)/rho1);
n = uw./ua;
a = 11e-3;
theta_droplet = [0:pi/1000:2*pi];

k = [1];
t = a/uw%(839-47)*0.05e-6

%% Ray tracing
alpha = [0:pi/1000:2*pi];
for n_level = 1:length(n)
for i = 1:length(alpha)
    if alpha(i) < asin(1/n(n_level))
        theta(i) = asin(n(n_level)*sin(alpha(i)));
        gamma(i) = 2*k*theta(i)-alpha(i)-(k-1)*pi;
        xm{n_level}(i) = (uw*t-n(n_level)*a*(1-cos(alpha(i)))-2*(k-1)*a*cos(theta(i)))*cos(gamma(i)-theta(i))-a*cos(gamma(i)-2*theta(i));
        ym{n_level}(i) = (uw*t-n(n_level)*a*(1-cos(alpha(i)))-2*(k-1)*a*cos(theta(i)))*sin(gamma(i)-theta(i))-a*sin(gamma(i)-2*theta(i));
    end
end
plot(xm{n_level},ym{n_level},'r'); hold on;
end

%% Plots
plot(a*cos(theta_droplet),a*sin(theta_droplet),'k'); hold on;
axis equal;
xlim([-a a])
ylim([-a a])

%% Fit circle

for i = 1:n_level
    clear R
    clear XC
    clear YC
    [R,XC,YC] = circfit(xm{i},ym{i});
    r(i) = R;
    xc(i) = XC;
    yc(i) = YC;
end

close all
plot(M,r./11e-3,'k')

% plot(xm{n_level},ym{n_level},'r'); hold on; axis equal;
% plot(XC+R*cos(alpha),YC+R*sin(alpha),'b')