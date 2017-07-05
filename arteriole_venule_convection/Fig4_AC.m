%%%
% Computes and plots results shown in panels Fig. 4, A and D, in Holter et al. (2017).
% Diffusion outside a cylinder, boundary is at cylinder wall, where
% concentration is set to C1. Initial condition is C1 at cylinder surface,
% C0 outside the cylinder, from the cylinder surface at r to infinity.
%
% Klas H. Pettersen.  
clear all; close all;
addpath('functions')
%addpath('functions/altmany-export_fig-f0af704')
fig_path = 'figs/';

%% Define parameters
D = [2.2e-6 1.96e-5 0.38e-6]; % cm^2/s [3kDa Texas red-dextran, potassium, dextran-70], dextrans from Sykova and Nicholson (2008), potassium from Halnes et al. (2013)
tortousity_vec = [2.04  1.6 2.13]; % same refs as for D
Dstar = D./tortousity_vec.^2; % cm^2/s [3kDa Texas red-dextran, potassium, large metabolite/dextran-70 (Asgari)]
tnum = [0 5 5*60 30*60]; % time points for concentration gradient plots
tmax = tnum(end); % maksimum sekunder
a = 23e-6*1e2/2; % m -> cm, radius of vessel. Asgari et al (2016)
nipol = 5000; % spatial resolution
r0 = a*1e-2;dmax = 200e-6; % r0 and dmax in meters


%% Concentration gradients at given times
figure(); % create figure
line_stl = {'--',':','-'};
xwidth=3;   % used in figure tp{kk} = plotting/saving (PaperSize and PaperPosition)
yheight=2;
hold on
for jj=1:3 % compute spatial consentration profile tnum time traces and three different diffusion constants
     tleg{jj} = ['D_{eff} = ' num2str(round(Dstar(jj),2,'significant')) ' \mu{}m^2/s']; % legend text
      
%    %% Parameters and numerical solution of differential equations
    xx = linspace(r0,dmax,nipol); % solve differential equation at these points

    xnum = unique(sort([xx 100e-6 300e-6 400e-6 800e-6 1500e-6 2500e-6 8000e-6])); % add some distant to increase accuracy
    u = pdex_test_v2(r0,dmax,Dstar(jj)*1e-4,xnum,tnum); % solve diff. eq., find concentration profile

    legend_ent = find_legend_text(tnum); % legend texts for all time traces
    speccols = [0 0 0]; speccols(jj) = 1;
 % Plot concentration profile 
    for ii=2:length(tnum)
        plot(xnum(1:(end-5))*1e3,100*u(ii,1:(end-5)),line_stl{ii-1},'LineWidth',1.5,'Color',speccols)
    end 
end
this_leg_num = legend(legend_ent,'FontSize',8,'Location','SouthEast');
xlabel('Position (mm)')
ylabel('Concentration (%)')
%print(gcf,'-dpdf',[fig_path 'diffusion_space_profile_cone_allDs.pdf'])
xlim([0 0.2]) % cm -> mm
ylim([0 100])
%legend(tleg,'FontSize',8,'Location','SouthEast');
set(gcf,'Units','Inches');
set(gcf, 'PaperSize', [xwidth yheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 xwidth yheight]);
print(gcf,'-dpdf',[fig_path 'diffusion_cone_space_profile.pdf'])


%% Concentration in a distance x as a function of time
figure()
hold on
for jj=1:3
    xnum2 = unique(sort([linspace(r0,200e-6,1000) 100e-6 300e-6 400e-6 800e-6 1500e-6 2500e-6]));
    tnum2 = [0 1/60 2/60 3/60 4/60 5/60 7/60 10/60 20/60 30/60 45/60 1 2 3 4 5 7.5 10 15 20 30 45 60]*60;
    u2 = pdex_test_v2(r0,dmax,Dstar(jj)*1e-4,xnum2,tnum2);
    cols = [0 0 0]; cols(jj) = 1;
    n100_2 = find(xnum2==100e-6); % m -> mm
    c100_2 = u2(:,n100_2);
    plot(tnum2(2:end)/60,100*c100_2(2:end),'LineWidth',1.5,'Color',cols)
end
xlim([0 60]) % cm -> mm
ylim([0 100])
legend(tleg,'FontSize',8,'Location','SouthEast');
xlabel('Time (min)')
ylabel('Concentration (%)')
set(gcf, 'PaperSize', [xwidth yheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 xwidth yheight]);
print(gcf,'-dpdf',[fig_path 'diffusion_cone100mum_v2.pdf'])