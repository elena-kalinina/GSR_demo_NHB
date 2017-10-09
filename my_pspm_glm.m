function scr_glm=my_pspm_glm(SubjID, sesN, cfg)

global settings;

if nargin < 3 
    cfg=[];
end


if ~isfield(cfg, 'PlotResults')
cfg.PlotResults=0; 

end

if ~isfield(cfg, 'DataPath')
cfg.DataPath='/home/elenka/Documents/MATLAB/MyGSRData/'; 

end


for i=1:length(sesN)
file_name=sprintf('%s%s%d.mat', cfg.DataPath, SubjID(1:2), sesN(i));
design_file=sprintf('%s%s_pspm_design_%d.mat', cfg.DataPath, SubjID(1:2), sesN(i));
fn{i}=file_name;
cfg.timing{i}=design_file;
end

%load the data
outfile=my_scr_import(fn, cfg);


%cfg.from=121;
%cfg.to=430;
%trim the data
outfile1=my_scr_trim(outfile, cfg);

%run glm on it
cfg.bf = settings.glm(1).cbf;
%Number of basis functions in the model apart from the canonical response
cfg.bf.args = 1;
[modelfile, scr_glm] = my_scr_glm(SubjID, outfile1, cfg);


if cfg.PlotResults
PlotN=[3, 5]; %1 design matrix 2 orthogonality 3 actual and predicted responses 4 regressors 5 reconstructed responses
fig = my_plot_scr_glm(modelfile, scr_glm, PlotN)
end

return;

