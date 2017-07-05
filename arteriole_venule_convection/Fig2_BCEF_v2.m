% Plots panels B, C, E, F in Fig. 2 in Holter et al. (2017).
% For B and C csv-files in folder 2c_slicedata_with_area is read.
% For E and F csv-files in folder 2e_slicedata_with_area is read.

% Klas H. Pettersen 03.07.2017, klas.pettersen@gmail.com

close all
clear all
% Load csv into vector of (dA, dv) tupples.
nbins = 1000; % bin resolution
tracer_depth = 100e-6; % check water arrival at this depth

geometry = '2c_slicedata_with_area';
%geometry = '2e_slicedata_with_area';

this_folder = ['./data/' geometry '/'];
figs_folder = ['./figs/'];

    
fnames = dir([this_folder '*.csv']);
for ii=1:length(fnames); % import all cross-sections (filenames)
    M=readtable(['./data/' geometry '/' fnames(ii).name]);
    element{ii}.dA = table2array(M(:,7)); % membrane area normal to x-component
    element{ii}.dv = table2array(M(:,4)); % x-component of velocity 
    legend_entries{ii} = fnames(ii).name(1:5);
    element{ii}.mean = sum(element{ii}.dA.*element{ii}.dv)/sum(element{ii}.dA);
    clear M
end
legend_entries{ii+1} = 'mean';

Nplanes = length(fnames);

above_x = 50;  
for ii=1:Nplanes
    disp(['Cross-section #' num2str(ii) ':'])
    disp(['Mean: ' num2str(element{ii}.mean) ' nm/s'])
    disp(['Net flow: ' num2str(sum(element{ii}.dA.*element{ii}.dv)*1e-9) ' mum^3/s'])
    [tmp ind] = sort(element{ii}.dv);
    element{ii}.dv = element{ii}.dv(ind);
    element{ii}.dA = element{ii}.dA(ind);
    maxv(ii) = max(element{ii}.dv);
    minv(ii) = min(element{ii}.dv);
    % Percentage of negative velocities:
    negvel = sum(element{ii}.dA(find(element{ii}.dv<0)))/sum(element{ii}.dA);
    disp([num2str(negvel*100) '% of the area contains negative velocities'])
    vel_above_x(ii) = sum(element{ii}.dA(find(element{ii}.dv>above_x)))/sum(element{ii}.dA);
    disp(['Net flow velocity above ' num2str(above_x) ' mum: ' num2str(vel_above_x(ii)*100) '%'])
    disp('--')
end
disp(['Mean net flow above ' num2str(above_x) ': ' num2str(mean(vel_above_x)*100) '%'])

min_value = min(0,min(minv)); % "global" minimum velocity
dv_bin_edges = min_value:(max(maxv)-min_value)/nbins:max(maxv); % same range for all traces
bin_size = dv_bin_edges(2)-dv_bin_edges(1);

A = zeros(length(dv_bin_edges)-1,Nplanes);
% Bin velocities:
for kk=1:Nplanes % all planes
    jj = 1;     % lower subscript (first bin)
    for ii=1:length(element{kk}.dv) % for all velocities
        while element{kk}.dv(ii)>dv_bin_edges(jj+1) % jump to next bin if next velocity is above lower bin size
                jj=jj+1; % move to next bin (next subscript)
        end
        A(jj,kk) = A(jj,kk)+element{kk}.dA(ii); % add area for velocity to its corresponding bin
    end
end


%% Panel B or E
fig_size = [3 2];
h1=figure('units','inches','position',[4 4 fig_size(1) fig_size(2)]);
set(gca,'fontsize',10)
x = dv_bin_edges(1:end-1)+bin_size/2; % position x-value at center of each bin
hold on
Atot = 0;
for ii=1:size(A,2) % all planes
    Atot = Atot+sum(A(:,ii)); % total area, all planes
end
Aav = Atot/size(A,2); % average cross sectional area
for ii=1:size(A,2)
    plot(x,A(:,ii)/(bin_size*Aav)); % area normalized to mean extracellular cross sectional area
end
plot(x,mean(A,2)/(bin_size*Aav),'LineWidth',1.5,'Color',[0 0 0])

xlabel('Velocity (nm/s)')
ylabel('Velocity density (norm.)')
xlim([-5 50])
ylim([0 0.25])
legend(legend_entries,'Location','NorthEast','fontsize',9)
set(h1,'Units','Inches');
pos = get(h1,'Position');
set(h1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1,[figs_folder geometry 'velocity_distribution.pdf'],'-dpdf','-r0')

disp(['Sum area: ' num2str(sum(mean(A,2)/(Aav)))]);

%figure() % Cumulative velocity
%TC = cumsum(mean(A,2)/(bin_size*Aav));
%plot(x,TC/TC(end))

%% Panel C or F
all_elements.dA = [];
all_elements.dv = [];
for ii=1:Nplanes
    all_elements.dA= [all_elements.dA; element{ii}.dA];
    all_elements.dv= [all_elements.dv; element{ii}.dv];
    element{ii}.times = tracer_depth./(element{ii}.dv*1e-9); % time of arrival at tracer_depth for elements in section ii
    temp = cumsum(fliplr(element{ii}.dA)); % cumulative sum of arrival times from low to high arrival times
    element{ii}.cum_sum = temp/sum(element{ii}.dA)*100; % fraction of area with arrived tracer at increasing times, normalized to percent
    clear temp
end

% Trace for mean:
[tmp ind] = sort(all_elements.dv);
all_elements.dv = all_elements.dv(ind);
all_elements.dA = all_elements.dA(ind);
all_elements.times = tracer_depth./(all_elements.dv*1e-9);
cum_sum_A = cumsum(fliplr(all_elements.dA));
cum_sum_A = cum_sum_A/sum(all_elements.dA)*100; % normalize to percent

h2=figure('units','inches','position',[5 5 fig_size(1) fig_size(2)]);
hold on
for ii=1:size(A,2) % all sections
    plot(flipud(element{ii}.times)/60,element{ii}.cum_sum);
end
plot(flipud(all_elements.times)/60,cum_sum_A,'LineWidth',1.5,'Color',[0 0 0]) % mean
set(gca,'fontsize',12)
xlim([0 1000])
ylim([0 100])
xlabel('Time (min)')
ylabel('Water (%)')
set(h2,'Units','Inches');
pos = get(h2,'Position');
set(h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2,[figs_folder geometry '_cumulative_velocity.pdf'],'-dpdf','-r0')
 
% Get average velocity distribution by taking the inverse of the avarage
% cumulative distribution?
