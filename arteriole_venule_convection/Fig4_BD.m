%%%
% This program reproduces the panels of Fig. 4 in Holter et al. (2017)
% for the planar geometry (panels B and D).
% Concentration profiles at different times and for solutes of different
% molecular weight (different diffusion coeffiscints) is plotted. 
%
% Klas H. Pettersen, 2017-03-30

clear all; close all;
addpath('functions')
fig_path = 'figs/';

% Parameters
my_pos = 0.01; % (cm) Plot concentration as a function of time at this position. 
D = [2.2e-6 1.96e-5 0.38e-6]; % cm^2/s [3kDa Texas red-dextran, potassium, large metabolite/dextran-70] % Sykova and Nicholson (first and last); Halnes et al (middle)
tortousity_vec = [2.04 1.6 2.13]; % Sykova and Nicholson (first and last); Halnes et al (middle)
Dstar = D./tortousity_vec.^2;
tmax = 60*60; % (s) maximal time in plot and in simulation
xmax = 0.02; % (cm) maximal spatial position in plot
x = 0:xmax/100:xmax; % position in cm
t_inst = [5 5*60 30*60];  % Plot the concentration profiles at these time instances.
c0=1; % initial concentration

%% Plot Fig. 4B
fig_size = [3 2];
figure('units','inches','position',[5 5 fig_size(1) fig_size(2)]);
hold on
tstyle = {'--',':','-'};
for jj=1:3
    for ii=1:length(t_inst)
        diff_erf{jj}(ii,:) = c0*erf(x/sqrt(4*Dstar(jj)*t_inst(ii)));
        ccol = [0 0 0];
        ccol(jj) = 1;
        plot(x*10,diff_erf{jj}(ii,:)*100,tstyle{ii},'LineWidth',1.5,'Color',ccol)
    end
end
xlabel('Position (mm)')
ylabel('Concentration (%)')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,[fig_path 'diffusion_space_profile_all.pdf'],'-dpdf','-r0')

%% Plot Fig. 4D
h2=figure('units','inches','position',[5 5 fig_size(1) fig_size(2)]);
hold on
t=0:tmax/200:tmax;
for ii = 1:length(Dstar)
    ccol = [0 0 0];
    ccol(ii) = 1;
    plot(t/60,100*(erf(my_pos./sqrt(4*Dstar(ii)*t))),'LineWidth',1.5,'Color',ccol)
    legend_d{ii} = ['D_{eff} = ' num2str(round(Dstar(ii)*1e8)) ' \mu{}m^2/s']; % cm^2/s -> mum^2/s
end
xlabel('Time (min)')
ylabel('Concentration (%)')
legend(legend_d,'Location','NorthEast','FontSize',8)
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2,[fig_path 'diffusion100mum.pdf'],'-dpdf','-r0')

%% Extra figures, with legends
dlabel = {'Texas Red','Potassium','70 kD dextran'};
for jj = 1:length(Dstar)
    h1(jj)=figure('units','inches','position',[5 5 fig_size(1) fig_size(2)]);
    tstyle = {'-','-','-'};
    hold on
    % Plot concentration profiles for time instances:
    for ii = 1:length(t_inst)
        ccol = [0 0 0];
        ccol(jj) = ii/length(t_inst);    % set color
        plot(x*10,diff_erf{jj}(ii,:)*100/c0,tstyle{jj},'LineWidth',1.5,'Color',ccol) % x*10: cm -> mm; diff_erf*100/c0: concentration-> per cent
        % Get right time unit in legend:
        if t_inst(ii)/60>1 % write minutes in legend
            ttext = num2str(round(t_inst(ii)/60*10)/10);
            utext = ' min';
        elseif round(t_inst(ii)*1000)==0 % use microsceonds
            ttext = num2str(round(t_inst(ii)*1000000));
            utext = ' \mu{}s';
        elseif round(t_inst(ii)*10)/10==0 % use milliseconds
            ttext = num2str(round(t_inst(ii)*1000));
            utext = ' ms';
        else % use seconds
           ttext = num2str(round(t(ii)*10)/10);
           utext = ' s';
        end 
        legend_ent{ii} = ['t = ' ttext utext];            
    end
    this_leg = legend(legend_ent,'FontSize',8);
    title(dlabel{jj})
end
