function legend_entry = find_legend_text(tnum)
t(1) = 0;
    for ii=2:length(tnum) % plot for all times tnum>0
        t(ii)=tnum(ii);
            if t(ii)/60>1 % write minutes in legend
                ttext = num2str(round(t(ii)/60*10)/10);
                utext = ' min';
            elseif round(t(ii)*1000)==0 % use microsceonds
                ttext = num2str(round(t(ii)*1000000));
                utext = ' \mu{}s';
            elseif round(t(ii)*10)/10==0 % use milliseconds
                ttext = num2str(round(t(ii)*1000));
                utext = ' ms';
            else % use seconds
               ttext = num2str(round(t(ii)*10)/10);
               utext = ' s';
            end 
 %           legend_ent{ii-1} = ['t = ' ttext utext];
            legend_entry{ii-1} = ['t = ' ttext utext];
       end
end

