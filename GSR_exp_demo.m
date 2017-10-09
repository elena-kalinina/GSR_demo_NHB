function GSR_exp_demo

%%EXPERIMENT SETTINGS
clear all

SubjID='vo';
ExpType='GSR_emo_RT';
numses=4;
error_all=zeros(1, numses);
means_all=zeros(1, numses);

addpath(genpath('PSPM_v3.0.2'));

cfg.experiment=1;

cfg.bl_dur=7;
cfg.async=1;
cfg.sr=10;
cfg.tau=cfg.sr;
cfg.tau_bar=1;
cfg.to_cut=60;
cfg.trial_length=10;
cfg.level=1.8; %1.2
cfg.factor=1;
cfg.baseline_mean=1;
cfg.rec_sr=49; %51.2;
cfg.Median=20;
cfg.settings.signal =1;
cfg.Window=5*round(cfg.sr)+1;
cfg.ForPSPM=1;


%% Import previous session data and protocols
sesN=3;
ExpType='GSR_emo_RT';
SubjID='vo';

GSR_import_RT_data(SubjID, ExpType,  sesN, cfg); 
get_protocol_data_test(SubjID, ExpType, sesN, cfg);

%% %% %% Running GLM

cfg.PlotResults=1; 

if sesN==1
filename=sprintf('%s_%s-%i.mat',SubjID, ExpType, sesN);
protocol=load(filename, 'onsets', 'control');
control=protocol.control';
onsets= protocol.onsets-cfg.to_cut;  
else
    filename=sprintf('%s_%s-%i.mat', SubjID, ExpType, sesN);
protocol=load(filename, 'onsets', 'control');
protocol_c=load(sprintf('%s_%s_c_%d.mat',  SubjID, ExpType, sesN));
onsets= protocol.onsets-cfg.to_cut;   
control=protocol_c.predicted_control;
end


if length(onsets)>20
    H=20;
else
    H=length(onsets);
end

H1=4; %6;
if cfg.experiment 
n_trial=H;
else
    n_trial=24; %or 22
end
n_trial
cfg.n_trial=n_trial;
cfg.n_tau=cfg.n_trial*cfg.trial_length;
sr=cfg.sr;
%sr=5;
tau=cfg.tau; %5; %secs
tau_bar=cfg.tau_bar; %sampling units

tr_length=cfg.trial_length; %secs
n_tau_bar=n_trial*tr_length*sr;
n_tau=n_trial*tr_length; 

cfg.from=cfg.to_cut; 
cfg.to=390; 

cfg.FeedbackFile=sprintf('%s_%s_c_%d.mat', SubjID, ExpType, sesN+1);
cfg.MatrixFile=sprintf('Gmatrix_%s_%d.mat', SubjID, sesN+1);

%[1,2];
[myglm, filter]=my_pspm_glm_filt(SubjID, sesN, cfg);

cfg.filter=filter;

observed_resp=myglm.Y;
%length(observed_resp)
predicted_resp_all=myglm.Yhat;
if myglm.recon(1)>0
cfg.param=myglm.recon(1);
else
cfg.param=abs(myglm.recon(1)); %1;
end
onsets=onsets*sr;  

onsets_neg=onsets(control(1:length(onsets))==1);
onsets_neutr=onsets(control(1:length(onsets))==0);

timeline1=zeros(n_tau_bar, 1); 

 for tm=1:length(timeline1)
     if ismember(tm, onsets_neg)
         timeline1(tm)=1;
     end
 end

observed_resp=observed_resp(onsets(1)-1:onsets(1)+n_tau_bar-1);  
length(observed_resp)

mycontrol= [observed_resp(1); timeline1(1:n_tau_bar-1)]; 
predicted_control=zeros(length(onsets), 1);
cfg.mymean=mean(observed_resp);
means_all(sesN)=mean(observed_resp);
cfg.mystd=std(observed_resp);
cfg.level=1.8; %1.2
cfg.factor=1; %0.8;

 if cfg.baseline_mean
     predicted_resp=mean(observed_resp(1:cfg.to_cut*sr))+(cfg.level+cfg.factor)*std(observed_resp(1:cfg.to_cut*sr));
 else
predicted_resp=cfg.mymean+(cfg.level+cfg.factor)*cfg.mystd; 
 end

[G_small, u_hat]=learning_LS_2G_new2_test(H1,n_tau, n_tau_bar, tr_length, tau, tau_bar, observed_resp, predicted_resp, mycontrol, cfg.param);

cfg.FeedbackFile=sprintf('%s_%s_c_%d.mat', SubjID, ExpType, sesN+1);
cfg.MatrixFile=sprintf('Gmatrix_%s_%d.mat', SubjID, sesN+1);
save(cfg.MatrixFile, 'G_small');
if cfg.experiment
    filename=cfg.FeedbackFile;
    save(filename, 'predicted_control');
end


%% Error analysis of the session

H1=4;
H2=3;
filename=sprintf('%s_%s-%i.mat',SubjID, ExpType, sesN);
protocol=load(filename, 'onsets');
real_onsets=protocol.onsets;
real_onsets=(real_onsets-protocol.onsets(1))*cfg.sr+1;
new_real_onsets=real_onsets(H1+1:H2:length(real_onsets));
if sesN==1
new_real_onsets=[1; new_real_onsets]
else
    new_real_onsets=[1 new_real_onsets]
end

ses_data=observed_resp((protocol.onsets(1)-cfg.to_cut):end);

length(ses_data); 

if cfg.baseline_mean
     target_resp=mean(observed_resp(1:cfg.to_cut*sr))+(cfg.level+cfg.factor)*std(observed_resp(1:cfg.to_cut*sr));
 else
target_resp=cfg.mymean+(cfg.level+cfg.factor)*cfg.mystd; 
 end
target_resp= repmat(target_resp, 1, length(ses_data)); 

error=(target_resp-(ses_data')).^2;  

figure;
plot(target_resp)
hold on
plot(ses_data); 
legend('target signal', 'observed signal')
title('Target signal (2.8 stds away from the baseline mean) vs observed signal')

error_all(sesN)=sum(error);

%% Real time session

SubjID='vo';
ExpType='GSR_emo_RT';
sesN=4; %2, 3
sr=10;
cfg.level=1.8; 
cfg.factor=1; 

cfg.OnsetFile=sprintf('%s_%s_o_%d.mat', SubjID, ExpType, sesN);

get_protocol_data_trd(SubjID, ExpType, sesN, cfg);
GSR_data_processing_RT_demo(SubjID, ExpType, sesN, cfg)


