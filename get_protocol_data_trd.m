function get_protocol_data_trd(SubjectID, ExpType, sessionN, cfg)

if nargin<4
    cfg=[];
end

if ~isfield(cfg, 'NofPicturesNeg') 
    cfg.NofPicturesNeg = 15;
end

if ~isfield(cfg, 'NofPicturesNeu') 
    cfg.NofPicturesNeu = 40;
end

if ~isfield(cfg, 'NofPicturesSes') 
    cfg.NofPicturesSes = 45;
end

if ~isfield(cfg, 'forPSPM') 
    cfg.forPSPM = 1;
end

if ~isfield(cfg, 'ProbControl') 
    cfg.ProbControl = 0.2; %or 0.3
end

if ~isfield(cfg, 'NumSessions') 
    cfg.NumSessions=1;
end

if ~isfield(cfg, 'nRepeat') 
    cfg.nRepeat=1;
end

if ~isfield(cfg, 'factorialStructure') 
    cfg.factorialStructure = [cfg.nRepeat, cfg.NofPicturesSes];
end

if ~isfield(cfg, 'factorNames') 
cfg.factorNames = [' NumberOfRepetitions ', 'NumberOfPictures'];
end

if ~isfield(cfg, 'targetPicture') 
    cfg.targetPicture = 1:cfg.NofPicturesSes;
end

if ~isfield(cfg, 'targetPictureNeu') 
    cfg.targetPictureNeu = 1:cfg.NofPicturesNeu*cfg.NumSessions;
end

if ~isfield(cfg, 'targetPictureNeg') 
    cfg.targetPictureNeg = cfg.NofPicturesNeu*cfg.NumSessions+1:(cfg.NofPicturesNeu+cfg.NofPicturesNeg)*cfg.NumSessions;
end

if ~isfield(cfg, 'leadInTimeSec') 
    cfg.leadInTimeSec=1;
end

if ~isfield(cfg, 'refreshHz') 
    cfg.refreshHz=60.3;
end

fname=sprintf('%s_%s-%i.trd', SubjectID, ExpType, sessionN);

fid = fopen(fname, 'rt');
counter = 0;
aline = fgetl(fid);
while ~feof(fid)
    aline = fgetl(fid);
    counter = counter + 1;
    numLine = str2num(aline);
    nElements = length(numLine);
    res(counter, 1:nElements) = numLine;
    
end
fclose(fid);
code = res(:, 1);
control=res(:, 3);
control=control(code>0);

%%%%%%%%%% ONSETS %%%%%%%%%%%%%%%%
onsets=res(:, 2);
onsets=onsets(code>0)+cfg.bl_dur+cfg.async;

filename=sprintf('%s_%s-%i.mat', SubjectID, ExpType, sessionN); 

save(filename, 'control', 'onsets')

if cfg.forPSPM
    if sessionN==1
        prepare_for_pspm(SubjectID,ExpType,  sessionN, cfg);
    else
        prepare_for_pspm_1(SubjectID,ExpType,  sessionN, cfg);
    end
end

return 