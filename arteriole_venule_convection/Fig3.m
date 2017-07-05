% Computes and plots velocities shown in Fig. 3 in Holter et al. (2017).
%
% Klas H. Pettersen 26.10.2016, klas.pettersen@gmail.com

close all
clear all
addpath('functions')
fig_path = 'figs/';
center_of_arterioles = [-280e-6/2 0]; % [x1 y1] in meter, keep x-value negative
center_of_venules = -center_of_arterioles; 
d = abs(center_of_arterioles(1)-center_of_venules(1));
diameters = 30e-6; % diameter of arteriole and venule [m],
% arteriole diameter is 23e-6 in Asgari et al (2016) and ref therein
% arteriole is 30e-6 in Jin et al (2016), vein 40e-6 in Jin et al (2016)
pressure_arteriole_mmHg = (d-diameters)*0.5e3; % pressure of 0.5 mmHg/mm per vessel = 1 mmHg/mm -> mmHg
pressure_arteriole = pressure_arteriole_mmHg*133.322; % mmHg -> Pa
permeability = [14.68 10.77]*(1e-9)^2; % [m^2], 14.68 nm^2 (highest) and 10.77 (lowest)
porosity = 0.2; % 20% extracellular volume
viscosity = 0.8e-3; % [Pa s]
K =  permeability/viscosity; % [m^2/Pa/s]

pressure_venule = -pressure_arteriole; % [Pa]
delta_p = pressure_arteriole-pressure_venule; % [Pa]
lambda = 2*pi*K*delta_p/log((1-2*d/diameters)^2); % [m^2/s]



%% Prepare plot
% Define line plot

% Define surface plot
%xmin = min([min(center_of_arterioles(:,1)) min(center_of_venules(:,1)) min(plot_line_at(:,1))]);
%xmax = max([max(center_of_arterioles(:,1)) max(center_of_venules(:,1)) max(plot_line_at(:,1))]);
%ymin = min([min(center_of_arterioles(:,2)) min(center_of_venules(:,2)) min(plot_line_at(:,2))]);
%ymax = max([max(center_of_arterioles(:,2)) max(center_of_venules(:,2)) max(plot_line_at(:,2))]);

xmin = center_of_arterioles(1);
xmax = center_of_venules(1);
ymin = xmin; ymax=xmax;
xrange = xmax-xmin;
yrange = ymax-ymin;
npoints = 400; % for line plot
plot_line_at = [linspace(0,0,npoints); linspace(2*ymin,2*ymax,npoints)]';

nx = 311; % x resolution
ny = 311; % y resolution
plot_xgrid = linspace(round((xmin-0.8*xrange)*1e5)/1e5,round((xmax+0.8*xrange)*1e5)/1e5,nx); % x-range for plot, rounded for every 10th of a micron
plot_ygrid = linspace(round((ymin-0.8*yrange)*1e5)/1e5,round((ymax+0.8*yrange)*1e5)/1e5,ny); % x-range for plot, rounded for every 10th of a micron
dx = plot_xgrid(2)-plot_xgrid(1);

%% Estimate line charge lambda based on pressure differences, distances and diameters
lambda_nu = lambda*1e-3*1e12*3600 % [m^2/s] -> pL/hour/micrometer: (1 m^2 = 1000 L/m = 1e-3 L/mum, pico = 1e-12, 3600 s/hour
disp(['A pressure difference of ' num2str(2*abs(pressure_arteriole_mmHg)) ' mmHg' ...
    ' and a distance of ' num2str(d*1e6) ' micrometers '])
disp(['between vessel centers gives an arterial efflux of ' num2str(round(lambda_nu(1),2,'significant')) ' pL/hour/micrometer vessel'])
disp(['for the largest permeability.' ...
    ' (' num2str(round(lambda_nu(2),2,'significant')) ' for the lowest permeability.)'])

%% Line plots
fsize = 6;

% Line plot along y-axis:
figure()
set(gca,'FontSize',fsize)
lcorr = 1.41*1.1;
tlwidth = 1;
xwidth=0.9*1.5/lcorr;   % used in figure (PaperSize, PaperPosition)
yheight=3/lcorr;
hold on;
this_y = plot_line_at(:,2);
plot(2*d*lambda(1)./(porosity*pi*(d^2+4*this_y.^2))*1e9,this_y*1e6,'LineWidth',tlwidth,'Color',[0 0 0])
hold on
plot(2*d*lambda(2)./(porosity*pi*(d^2+4*this_y.^2))*1e9,this_y*1e6,'LineWidth',tlwidth,'Color',[0.6 0.6 0.6])
xlabel('Velocity [nm/s]')
ylabel('y [\mu{}m]')
set(gca,'YTick',[-300 -200 -100 0 100 200 300])
ylim([min(plot_ygrid)*1e6 max(plot_ygrid)*1e6])
set(gcf, 'PaperSize', [xwidth yheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 xwidth yheight]);
%print(gcf,'-dpdf',['poff' num2str(teller) 'temp.pdf']);
leg_array = {['\kappa = ' num2str(round(permeability(1)*1e18,2,'significant')) ' nm^2'],['\kappa = ' num2str(round(permeability(2)*1e18,2,'significant')) ' nm^2']};
%legend(leg_array)
print(gcf,'-dpdf',[fig_path 'alongyaxis.pdf']);

% Line plot along x-axis
figure()
set(gca,'FontSize',fsize)
xwidth=3.22/lcorr;   % used in figure (PaperSize, PaperPosition)
yheight=1.1/lcorr;
hold on;
this_x_all = plot_line_at(:,2);
lcols = [1 0 0;
         1 .6 .6];
for ii=1:2 % both permeabilities
this_x = sort(this_x_all(find(this_x_all<-d/2-diameters/2)));
h(ii)=plot(this_x*1e6,abs(2*d*lambda(ii)./((d^2*pi-4*pi*this_x.^2)*porosity)*1e9),'LineWidth',tlwidth,'Color',lcols(ii,:)); % vx
this_x = sort(intersect(this_x_all(find(this_x_all>-d/2+diameters/2)),this_x_all(find(this_x_all<d/2-diameters/2))));
plot(this_x*1e6,abs(2*d*lambda(ii)./((d^2*pi-4*pi*this_x.^2)*porosity)*1e9),'LineWidth',tlwidth,'Color',lcols(ii,:))
this_x = sort(this_x_all(find(this_x_all>d/2+diameters/2)));
plot(this_x*1e6,abs(2*d*lambda(ii)./(porosity*(d^2*pi-4*pi*this_x.^2))*1e9),'LineWidth',tlwidth,'Color',lcols(ii,:))
end
xlim([min(plot_xgrid)*1e6 max(plot_xgrid)*1e6])
ylim([0 40])
ylabel('Velocity [nm/s]')
xlabel('y [\mu{}m]')
set(gcf, 'PaperSize', [xwidth yheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 xwidth yheight]);
print(gcf,'-dpdf',[fig_path 'alongxaxis.pdf']);



%% Compute and plot velocity at grid
vx = zeros(length(plot_ygrid),length(plot_xgrid));
vy = zeros(length(plot_ygrid),length(plot_xgrid));
x = plot_xgrid;
for jj=1:length(plot_xgrid)
    for ii=1:length(plot_ygrid)
        x = plot_xgrid(jj);
        y = plot_ygrid(ii);
        if sqrt((x-d/2)^2+y^2)<(diameters/2) || sqrt((x+d/2)^2+y^2)<(diameters/2)  % inside arteriole or venule
            vx(ii,jj) = 0;
            vy(ii,jj) = 0;
        else
            vx(ii,jj) = (2*d*(d^2-4*x.^2+4*y.^2)*lambda(1))/(pi*porosity*(d^4+8*d^2*(-x.^2+y.^2)+16*(x.^2+y.^2).^2));
            vy(ii,jj) =  -((16*d*x*y*lambda(1))/(pi*porosity*(d^4+8*d^2*(-x.^2+y.^2)+16*(x.^2+y.^2).^2)));
        end
    end
end

% Plot
figure()
set(gca,'FontSize',fsize)
xwidth=3*1.5/lcorr;   % used in figure (PaperSize, PaperPosition)
yheight=2*1.5/lcorr;
clims=[0 max(max(sqrt(vx.^2+vy.^2)*1e9))];
imagesc(plot_xgrid*1e6,plot_ygrid*1e6,sqrt(vx.^2+vy.^2)*1e9,clims)
set(gca,'FontSize',fsize)
xlabel('x [\mu m]')
ylabel('y [\mu m]')
c=colorbar('FontSize',fsize);
ylabel(c,'Velocity [nm/s]')
hold on
axis equal
% Add grey and black line for line plots, and circles for arterioles and
% venules:
plot(plot_line_at(:,1)*1e6,plot_line_at(:,2)*1e6,'LineWidth',1.5,'Color',[0 0 0])
this_x_all = plot_line_at(:,2);
this_x = sort(this_x_all(find(this_x_all<-d/2-diameters/2)));
plot(this_x*1e6,zeros(1,length(this_x)),'LineWidth',1.5,'Color',[1 0 0])
this_x = sort(intersect(this_x_all(find(this_x_all>-d/2+diameters/2)),this_x_all(find(this_x_all<d/2-diameters/2))));
plot(this_x*1e6,zeros(1,length(this_x)),'LineWidth',1.5,'Color',[1 0 0])
this_x = sort(this_x_all(find(this_x_all>d/2+diameters/2)));
plot(this_x*1e6,zeros(1,length(this_x)),'LineWidth',1.5,'Color',[1 0 0])
filledCircle(center_of_arterioles*1e6,(diameters/2+dx/2)*1e6,200,'r');
filledCircle(center_of_venules*1e6,(diameters/2+dx/2)*1e6,200,'b');
xlim([plot_xgrid(1)*1e6 plot_xgrid(end)*1e6])
ylim([plot_ygrid(1)*1e6 plot_ygrid(end)*1e6])
set(gcf, 'PaperSize', [xwidth yheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 xwidth yheight]);
print(gcf,'-dpdf',[fig_path 'velocity_color.pdf']);