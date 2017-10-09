function ampl = calculate_amplitude(data, onsets, trial_len)

counter=0;
for ons = onsets
    counter=counter+1
    trial=data(ons:ons+trial_len);
    length(trial);
    max_ind=find(trial==max(trial));
    min_ind=find(trial==min(trial));
    
    if max_ind > min_ind
    
    trial_min = min(data(1:max_ind));
%     min_ind=find(trial==trial_min)
%     data(max_ind)
    ampl_trial(counter)=abs(data(max_ind)-trial_min);
    
    else
        
        trial_max=max(data(min_ind:end));
        ampl_trial(counter)=abs(trial_max-data(min_ind));
    end
    
    
end

ampl=mean(ampl_trial);

return 