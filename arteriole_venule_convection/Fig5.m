%%%
% This program plots Fig. 5 in Holter et al. (2017)
% Peclet numbers. 
%
% Klas H. Pettersen, 2017-03-30

close all;
clear all;
fig_path = 'figs/';

arteriole_diameter = 36e-6;
venule_diameter = 47.3e-6;
L = 280e-6-(arteriole_diameter+venule_diameter)/2;
permeability = 14.68*(1e-9)^2; viscosity = 0.8e-3;Vex = 0.2; dpdx = 133*1e3*[0.5 1 2];
u=permeability/viscosity/Vex*dpdx;
p2 = 10.77*(1e-9)^2; u2=p2/viscosity/Vex*dpdx;
disp(['Velocity 1: ' num2str(u*1e9) ' nm/s (for permeability of ' num2str(permeability) 'nm^2).'])
disp(['Velocity 2: ' num2str(u2*1e9) ' nm/s (for permeability of ' num2str(p2) 'nm^2).'])

solute = {'Dex70'   ,'Ovalbumin'    ,'Dex3',    'Sucrose',  'Na^+',  'Ca^{2+}'    ,'Dopamine' ,'Penicillin'   , 'A{\beta} monomer'};
D_exp= [0.38e-6    1.0e-6             2.2e-6       7.0e-6      20e-6   9.90e-6     6.9e-6      3.52e-6         1.8e-6]*1e-4; 
lambda  = [2.13      2.5             2.04        1.6         1.82    1.82        1.8         1.62                1.7];
D_eff = D_exp./lambda.^2;
T = table(solute',D_eff');
T.Properties.VariableNames = {'Solute' 'Diff_coeff'};
disp(T)
D = linspace(min(D_eff),max(D_eff),200);
xwidth=6;   % used in figure tp{kk} = plotting/saving (PaperSize and PaperPosition)
yheight=4;

figure()
for ii=1:3
    Pe(ii,:) = L*u(ii)./D_eff;
end
loglog(D_eff,Pe,'LineWidth',2)
hold on
%loglog(D,0.05*ones(1,length(D)),'--','Color',[0 0 0])
loglog(D_eff,Pe,'.')
for ii=1:length(D_eff)
    txt = ['\leftarrow ' solute{ii}];
    text(D_eff(ii),Pe(2,ii),txt);
end
xlabel('D^* (m^2/s)')
ylabel('Peclet number')
xlim([min(D_eff)/2 2*max(D_eff)]) % cm -> mm
%ylim([0 100])
set(gcf, 'PaperSize', [xwidth yheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 xwidth yheight]);
print(gcf,'-dpdf',[fig_path 'Peclet_numbers.pdf'])