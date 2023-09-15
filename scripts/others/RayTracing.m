%-------------------------------------------------------------------------------------------------------------------------------
%                                                           OBJECTIVE                                                           
%-------------------------------------------------------------------------------------------------------------------------------
% This code simulates the position of waves within a droplet. The latter being impacted by a shock wave.
% It corresponds to the analytical wave description described in:
% L. Biasiori-Poulanges, K. Schmidmayer. A phenomenological analysis of droplet shock-induced cavitation using a multiphase
% modelling approach. Physics of Fluids, Vol. 35, 013312, 2023.

%-------------------------------------------------------------------------------------------------------------------------------
%                                                             CODE                                                           
%-------------------------------------------------------------------------------------------------------------------------------

clear all; close all; clc;

%% Initial conditions
Ms = 2.4;

a = 11.0e-3;
b = 11.0e-3;
epsilon=(b/a)^2-1;
alpha=linspace(-pi/2,pi/2,1000);

c = 345;
uw = sqrt(4.4*(6e8+1.01325e5)/998);
ua = Ms*c;
n = uw/ua;

t_impact = (22e-3-20e-3)/ua;
Delta_frame = 0.5e-6;
t_frame5 = 5*Delta_frame;

index = 0;
for i = 28:12:75
    index = index+1;
    time(index) = (t_frame5-t_impact)+i*0.5e-6;
end
disp(time);

for j = 1:length(time)
    t=time(j);

    %% Wavefronts
    for i=1:length(alpha)
        
        delta(i) = -alpha(i)+asin(n*sin(alpha(i)));
        
        xb(i) = -a^2/sqrt(a^2+b^2*tan(alpha(i))^2);
        yb(i) = tan(alpha(i))*b^2/sqrt(a^2+b^2*tan(alpha(i))^2);
        
        xc(i) = xb(i)+(uw*t-n*(a-a^2/sqrt(a^2+tan(alpha(i))^2*b^2)))*cos(delta(i));
        yc(i) = yb(i)+(uw*t-n*(a-a^2/sqrt(a^2+tan(alpha(i))^2*b^2)))*sin(delta(i));
        
        xd(i) = (a^2*xb(i)*sin(delta(i))^2-b^2*xb(i)*cos(delta(i))^2-...
            a^2*yb(i)*sin(2*delta(i)))/(b^2*cos(delta(i))^2+a^2*sin(delta(i))^2);
        yd(i) = (b^2*yb(i)*cos(delta(i))^2-a^2*yb(i)*sin(delta(i))^2-...
            b^2*xb(i)*sin(2*delta(i)))/(b^2*cos(delta(i))^2+a^2*sin(delta(i))^2);
        
        ErefEx(i) = ((a^4*yd(i)^2-b^4*xd(i)^2)*cos(delta(i))-...
            2*a^2*b^2*xd(i)*yd(i)*sin(delta(i)))/(a^4*yd(i)^2+b^4*xd(i)^2);
        
        ErefEy(i) = ((b^4*xd(i)^2-a^4*yd(i)^2)*sin(delta(i))-...
            2*a^2*b^2*xd(i)*yd(i)*cos(delta(i)))/(a^4*yd(i)^2+b^4*xd(i)^2);
        
        l(i) = uw*t-n*(a+xb(i))+2*(b^2*xb(i)*cos(delta(i))+a^2*yb(i)*sin(delta(i)))/(b^2*cos(delta(i))^2+a^2*sin(delta(i))^2);
        
        xe(i) = xd(i)+l(i)*ErefEx(i);
        ye(i) = yd(i)+l(i)*ErefEy(i);
      
    end
    
    %% Droplet
    C = [0 0];  % center 
    th = linspace(0,2*pi) ; 
    x_droplet = C(1)+a*cos(th) ; 
    y_droplet = C(2)+b*sin(th) ;

    %% Graph Setting
    xcolor = linspace(-2*pi,2*pi);
    ycolor = sin(xcolor);

    str_green = '#008080';
    str_red = '#C60C1F';
    str_blue = '#050572';
    str_gray = '#3E3E40';

    color_green = sscanf(str_green(2:end),'%2x%2x%2x',[1 3])/255;
    color_red = sscanf(str_red(2:end),'%2x%2x%2x',[1 3])/255;
    color_blue = sscanf(str_blue(2:end),'%2x%2x%2x',[1 3])/255;
    color_gray = sscanf(str_gray(2:end),'%2x%2x%2x',[1 3])/255;

    %% Plots
    plot(xc,yc,'color',color_red,'linewidth',0.5); hold on; 
    plot(xe,ye,'color',color_green,'linewidth',0.5);

end

%Droplet
plot(x_droplet,y_droplet,'color',color_gray,'linewidth',2); 
axis equal


xlim([-a,a]);
ylim([-b,b]);
% xlabel('$x$~[m]', 'Interpreter', 'Latex', 'FontSize', 16)
% ylabel('$y$~[m]', 'Interpreter', 'Latex', 'FontSize', 16)

% set(0, 'defaultAxesTickLabelInterpreter','latex');
% set(0, 'defaultLegendInterpreter','latex');

aa = get(gca,'XTickLabel');
set(gca,'fontsize',16);

pause;