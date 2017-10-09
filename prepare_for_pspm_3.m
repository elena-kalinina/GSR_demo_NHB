function prepare_for_pspm_3(SubjID, ExpType, sesN, cfg)

names={'neg', 'neutr'};
load(sprintf('%s_%s-%d.mat', SubjID,ExpType, sesN), 'control');
load(sprintf('%s_%s-%d.mat', SubjID,ExpType, sesN), 'onsets');
neg_onsets=onsets(control(1:length(onsets))==1)-cfg.to_cut; %67; %onsets(1)+2; %+1; 
neutr_onsets=onsets(control(1:length(onsets))==0)-cfg.to_cut; %67; %onsets(1)+2; %+1; 
onsets=[];
onsets={neg_onsets, neutr_onsets}; 
%durations={3}; 
filename_save=sprintf('%s_pspm_design_%d.mat', SubjID, sesN);
save(filename_save, 'names', 'onsets')