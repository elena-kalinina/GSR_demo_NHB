function [modelfile, filter, glm] = my_scr_glm_filt(SubjID, datafile, cfg)
% scr_glm specifies a within subject general linear convolution model of 
% predicted signals and calculates amplitude estimates for these responses 
%
% FORMAT:
% glm = scr_glm(model, options)
%
% MODEL with required fields
% model.modelfile:  a file name for the model output
% model.datafile:   a file name (single session) OR
%                   a cell array of file names
% model.timing:     a multiple condition file name (single session) OR
%                   a cell array of multiple condition file names OR
%                   a struct (single session) with fields .names, .onsets, 
%                       and (optional) .durations and .pmod  OR
%                   a cell array of struct
% model.timeunits:  one of 'seconds', 'samples', 'markers'
% 
% optional fields
% model.modality:   currently 'scr', default 'scr'
% model.bf:         basis function/basis set; modality specific default
%                   with subfields .fhandle (function handle or string) and
%                   .args (arguments, first argument sampling interval will
%                   be added by scr_glm)
% model.channel:    channel number; default: first channel of the specified modality
% model.norm:       normalise data; default 0
% model.filter:     filter settings; modality specific default
% model.missing:    allows to specify missing (e. g. artefact) epochs in
%                   the data file. See scr_get_timing for epoch definition; 
%                   specify a cell array for multiple input files. This
%                   must always be specified in SECONDS.
%                   Default: no missing values
% model.nuisance:   allows to specify nuisance regressors. Must be a file
%                   name; the file is either a .txt file containing the
%                   regressors in columns, or a .mat file containing the
%                   regressors in a matrix variable called R. There must be 
%                   as many values for each column of R as there are data 
%                   values. SCRalyze will call these regressors R1, R2, ...

%%MODEL with required fields
%model.modelfile:  a file name for the model output
%model.datafile:   a file name (single session) OR
%                  a cell array of file names
%model.timing:     a multiple condition file name (single session) OR
%                   a cell array of multiple condition file names OR
%                   a struct (single session) with fields .names, .onsets, 
%                       and (optional) .durations and .pmod  OR
%                   a cell array of struct

if nargin < 3
    cfg= [];
end

if ~isfield(cfg, 'modelfile')
cfg.modelfile= sprintf('%s_model.mat', SubjID); % 'samples', 'markers'
end

if ~isfield(cfg, 'datafile')
cfg.datafile= datafile; %{sprintf('%stscr_%s%d.mat', cfg.DataPath, SubjID, 1), sprintf('%stscr_%s%d.mat', cfg.DataPath, SubjID, 2) }; % 'samples', 'markers'
end

if ~isfield(cfg, 'timing')
%cfg.timing= {sprintf('%s%s_pspm_design_%d.mat', cfg.DataPath, SubjID, 1), sprintf('%s%s_pspm_design_%d.mat', cfg.DataPath, SubjID, 2)}; 

cfg.timing= {sprintf('%s_pspm_design_%d.mat', SubjID, 1), sprintf('%s_pspm_design_%d.mat', SubjID, 2)};
% 'samples', 'markers'
end


if ~isfield(cfg, 'timeunits')
cfg.timeunits= 'seconds'; % 'samples', 'markers'
end


% OPTIONS (optional argument)
% options.overwrite: overwrite existing model output; default 0
% options.marker_chan_num: marker channel number; default first marker
%                          channel
%
% TIMING - multiple condition file(s) or struct variable(s):
% The structure is equivalent to SPM2/5/8/12 (www.fil.ion.ucl.ac.uk/spm), 
% such that SPM files can be used.
% The file contains the following variables:
% - names: a cell array of string for the names of the experimental
%   conditions
% - onsets: a cell array of number vectors for the onsets of events for
%   each experimental condition, expressed in seconds, marker numbers, or
%   samples, as specified in timeunits
% - durations (optional, default 0): a cell array of vectors for the 
%   duration of each event. You need to use 'seconds' or 'samples' as time 
%   units
% - pmod: this is used to specify regressors that specify how responses in
%   an experimental condition depend on a parameter to model the effect
%   e.g. of habituation, reaction times, or stimulus ratings.
%   pmod is a struct array corresponding to names and onsets and containing
%   the fields
%   - name: cell array of names for each parametric modulator for this
%       condition
%   - param: cell array of vectors for each parameter for this condition,
%       containing as many numbers as there are onsets
%   - poly (optional, default 1): specifies the polynomial degree
%
% e.g. produce a simple multiple condition file by typing
%  names = {'condition a', 'condition b'};
%  onsets = {[1 2 3], [4 5 6]};
%  save('testfile', 'names', 'onsets');
%
%
% RETURNS a structure 'glm' which is also written to file
%
% -------------------------------------------------------------------------
% REFERENCES:
%
% (1) GLM for SCR:
% Bach DR, Flandin G, Friston KJ, Dolan RJ (2009). Time-series analysis for 
% rapid event-related skin conductance responses. Journal of Neuroscience
% Methods, 184, 224-234.
%
% (2) Canonical response function, and GLM assumptions for SCR:
% Bach DR, Flandin G, Friston KJ, Dolan RJ (2010). Modelling event-related 
% skin conductance responses. International Journal of Psychophysiology,
% 75, 349-356.
% 
% (3) Fine-tuning of filters and response functions:
% Bach DR, Friston KJ, Dolan RJ (2013). An improved algorithm for
% model-based analysis of evoked skin conductance responses. Biological
% Psychology, 94, 490-497.
% 
%__________________________________________________________________________
% PsPM 3.0
% (C) 2008-2015 Dominik R Bach (Wellcome Trust Centre for Neuroimaging)

% $Id: scr_glm.m 701 2015-01-22 14:36:13Z tmoser $
% $Rev: 701 $

% initialise & user output
% -------------------------------------------------------------------------
global settings;
if isempty(settings), scr_init; end;
glm = struct([]); % output model structure
tmp = struct([]); % temporary model structure

% check input arguments & set defaults
% -------------------------------------------------------------------------

% check missing input --
if nargin<1
    errmsg=sprintf('Nothing to do.'); warning('ID:invalid_input', errmsg); return;
elseif nargin<5
    options = struct();
end;

fprintf('Computing GLM: %s ...\n', cfg.modelfile);

if ~isfield(cfg, 'datafile')
    warning('ID:invalid_input', 'No input data file specified.'); return;
elseif ~isfield(cfg, 'modelfile')
    warning('ID:invalid_input', 'No output model file specified.'); return;
elseif ~isfield(cfg, 'timing')
    warning('ID:invalid_input', 'No event onsets specified.'); return;
elseif ~isfield(cfg, 'timeunits')
    warning('ID:invalid_input', 'No timeunits specified.'); return;      
end;

% check faulty input --
if ~ischar(cfg.datafile) && ~iscell(cfg.datafile)
    warning('ID:invalid_input', 'Input data must be a cell or string.'); return;
elseif ~ischar(cfg.modelfile)
    warning('ID:invalid_input', 'Output model must be a string.'); return;
elseif ~ischar(cfg.timing) && ~iscell(cfg.timing) && ~isstruct(cfg.timing)
    warning('ID:invalid_input', 'Event onsets must be a string, cell, or struct.'); return;
elseif ~ischar(cfg.timeunits) || ~ismember(cfg.timeunits, {'seconds', 'markers', 'samples'})
    warning('ID:invalid_input', 'Timeunits (%s) not recognised; only ''seconds'', ''markers'' and ''samples'' are supported', cfg.timeunits); return;
end;

% get further input or set defaults --
% check modality --
if ~isfield(cfg, 'modality')
    cfg.modality = 'scr';
elseif ~ismember(cfg.modality, {settings.glm.modality})
    warning('ID:invalid_input', 'Unknown modality %s.', cfg.modality); return;
end;
modno = find(strcmpi(cfg.modality, {settings.glm.modality}));

% check data channel --
if ~isfield(cfg, 'channel')
    cfg.channel = cfg.modality; % or 1; this returns the first channel of this type
elseif ~isnumeric(cfg.channel)
    warning('ID:invalid_input', 'Channel number must be numeric.'); return;
end;

% check normalisation --
if ~isfield(cfg, 'norm')
    cfg.norm = 0;
elseif ~ismember(cfg.norm, [0, 1])
    warning('ID:invalid_input', 'Normalisation must be specified as 0 or 1.'); return; 
end;

% check filter --
if ~isfield(cfg, 'filter')
    cfg.filter = settings.glm(modno).filter;
    filter=cfg.filter;
elseif ~isfield(cfg.filter, 'down') || ~isnumeric(cfg.filter.down) 
    % tested because the field is used before the call of scr_prepdata (everything else is tested there)
    warning('ID:invalid_input', 'Filter structure needs a numeric ''down'' field.'); return;
else
    filter=cfg.filter;
end;

% check & get basis functions --
basepath = [];
if ~isfield(cfg, 'bf')
    cfg.bf = settings.glm(modno).cbf;
else
    if ~isfield(cfg.bf, 'fhandle')
        warning('No basis function given.');
    elseif ischar(cfg.bf.fhandle)
        [basepath, basefn, baseext] = fileparts(cfg.bf.fhandle);
        cfg.bf.fhandle = str2func(basefn);
    elseif ~isa(cfg.bf.fhandle, 'function_handle')
        warning('Basis function must be a string or function handle.');
    end;
    if ~isfield(cfg.bf, 'args')
        cfg.bf.args = [];
    elseif ~isnumeric(cfg.bf.args)
        warning('Basis function arguments must be numeric.');
    end;
end;
if ~isempty(basepath), addpath(basepath); end; 
try
    cfg.bf.X = feval(cfg.bf.fhandle, [1/cfg.filter.down; cfg.bf.args(:)]);
catch
    warning('ID:invalid_fhandle', 'Specified basis function %s doesn''t exist or is faulty', func2str(cfg.bf.fhandle)); return;
end;

% remove path & clear local variables --
if ~isempty(basepath), rmpath(basepath); end;
clear basepath basefn baseext
% check options --
if ~isfield(options, 'overwrite')
    options.overwrite = 0;
elseif ~ismember(options.overwrite, [0, 1])
    options.overwrite = 0;
end; 
if ~isfield(options, 'marker_chan_num')
    options.marker_chan_num = 'marker';
elseif ~(isnumeric(options.marker_chan_num) && numel(options.marker_chan_num)==1)
    options.marker_chan_num = 'marker';
end; 

% check files --
if exist(cfg.modelfile, 'file') && ~(isfield(options, 'overwrite') && options.overwrite == 1)
    overwrite=menu(sprintf('Model file (%s) already exists. Overwrite?', cfg.modelfile), 'yes', 'no');
    if overwrite == 2, return, end;
    options.overwrite = 1;
end;

if ischar(cfg.datafile)
    cfg.datafile={cfg.datafile};
end;
if ischar(cfg.timing) || isstruct(cfg.timing)
    cfg.timing = {cfg.timing};
end;

if numel(cfg.datafile) ~= numel(cfg.timing)
    warning('ID:number_of_elements_dont_match', 'Session numbers of data files and event definitions do not match.'); return;
end;

% check & get data --
fprintf('Getting data ...');
nFile = numel(cfg.datafile);
for iFile = 1:nFile
    [sts, infos, data] =  my_scr_load_data(cfg.datafile{iFile}, cfg.channel); %scr_load_data(cfg.datafile{iFile}, cfg.channel);
    if sts < 1, return; end;
    y{iFile} = data{1}.data(:);
    sr(iFile) = data{1}.header.sr;
    fprintf('.');
    if any(strcmp(cfg.timeunits, {'marker', 'markers'}))
        [sts, infos, data] = scr_load_data(cfg.datafile{iFile}, options.marker_chan_num);
        if sts < 1, return; end;
        events{iFile} = data{1}.data * data{1}.header.sr;
    end;
end;
if nFile > 1 && any(diff(sr) > 0)
    fprintf('\nSample rate differs between sessions.\n')
else
    fprintf('\n');
end;

% check regressor files --
[sts, multi] = scr_get_timing('onsets', cfg.timing, cfg.timeunits);
if sts < 0, warning('Invalid multiple condition file'); return; end;

% check & get missing values --
if ~isfield(cfg, 'missing')
    missing = cell(nFile, 1);
else
    if ischar(cfg.missing) || isnumeric(cfg.missing)
        cfg.missing = {cfg.missing};
    elseif ~iscell(cfg.missing)
        warning('ID:invalid_input', 'Missing values must be a filename, matrix, or cell array of these.'); return;
    end;
    if numel(cfg.missing) ~= nFile
        warning('ID:number_of_elements_dont_match', 'Same number of data files and missing value definitions is needed.'); return;
    end;
    for iSn = 1:nFile
        if isempty(cfg.missing{iSn})
            sts = 1; missing{iSn} = [];
        else
            [sts, missing{iSn}] = scr_get_timing('epochs', cfg.missing{iSn}, 'seconds');
        end;
        if sts == -1, return; end;
    end;
end;

% check and get nuisance regressors
if ~isfield(cfg, 'nuisance')
    cfg.nuisance = cell(nFile, 1);
    for iSn = 1:nFile
        R{iSn} = [];
    end
    nR = 0;
else
    if ischar(cfg.nuisance)
        cfg.nuisance = {cfg.nuisance};
    elseif ~iscell(cfg.nuisance)
        warning('ID:invalid_input', 'Nuisance regressors must be specified as char or cell of file names.'); return;
    end;
    if numel(cfg.nuisance) ~= nFile
        warning('ID:number_of_elements_dont_match', 'Same number of data files and nuisance regressor files is needed.'); return;
    end;

    for iSn = 1:nFile
        if isempty(cfg.nuisance{iSn})
            R{iSn} = [];
        else
            try
                indata = load(cfg.nuisance{iSn});
                if isstruct(indata)
                    R{iSn} = indata.R;
                else
                    R{iSn} = indata;
                end;
            catch
                warning('ID:invalid_file_type', 'Unacceptable file format or non-existing file for nuisance file in session %01.0f', iSn); return;
            end;
            if size(R{iSn}, 1) ~= numel(y{iSn})
                warning('ID:number_of_elements_dont_match', 'Nuisance regressors for session %01.0f must have same number of data points as observed data.', iSn); return;
            end;
        end;
        if iSn == 1
            nR = size(R{iSn}, 2);
        elseif size(R{iSn}, 2) ~= nR
            warning('ID:number_of_elements_dont_match', 'Nuisance regressors for all sessions must have the same number of columns'); return;
        end;
    end;
end


fprintf('Preparing & inverting model ... ');

% collect output model information --
glm(1).glmfile    = cfg.modelfile; % this field will be removed in the future so don't use any more
glm.modelfile     = cfg.modelfile;
glm.input         = cfg; %model;
glm.input.options = options;
glm.bf            = cfg.bf;
glm.bf.bfno       = size(glm.bf.X, 2);

% clear local variables --
clear sts iFile modno


% prepare data & regressors
%-------------------------------------------------------------------------

Y=[]; M=[]; tmp=struct([]);
for iSn = 1:nFile
    
    % prepare (filter & downsample) data
    cfg.filter.sr = sr(iSn);
    [sts, newy, newsr] = scr_prepdata(y{iSn}, cfg.filter);
    if sts ~= 1, return; end;
    
    % concatenate data 
    Y=[Y; newy(:)];
    
    % get duration of single sessions
    tmp(1).snduration(iSn) = numel(newy);
    
    % process missing values
    newmissing = zeros(size(newy(:)));
    if ~isempty(missing{iSn})
        missingtimes = missing{iSn} * newsr;
        for iMs = 1:size(missingtimes, 1)
            newmissing(missingtimes(iMs, 1):missingtimes(iMs, 2)) = 1;
        end;
    end;
    M = [M; newmissing];    
       
    % convert regressor information to samples
    for n = 1:numel(multi(1).names)
        % convert onsets to samples
        switch cfg.timeunits
            case 'samples'
                newonsets    = round(multi(iSn).onsets{n} * newsr/sr(iSn));
                newdurations = round(multi(iSn).durations{n} * newsr/sr(iSn));
            case 'seconds'
                newonsets    = round(multi(iSn).onsets{n} * newsr);
                newdurations = round(multi(iSn).durations{n} * newsr);
            case 'markers'
                try
                    newonsets = round(events{iSn}(multi(iSn).onsets{n}) * newsr); % markers are timestamps in seconds
                catch
                    warning('\nSome events in condition %01.0f were not found in the data file %s', n, cfg.datafile{iSn}); return;
                end;
                newdurations = multi(iSn).durations{n};
        end;
        % get the first multiple condition definition --
        if iSn == 1
            names{n} = multi(1).names{n};
            onsets{n} = [];
            durations{n} = [];
            if isfield(multi, 'pmod') && (numel(multi(1).pmod) >= n)
                for p = 1:numel(multi(1).pmod(n).param)
                    pmod(n).param{p} = [];
                end
                pmod(n).name = multi(1).pmod(n).name;
            end;
        % or shift multiple condition definition --
        else
            newonsets = newonsets + sum(tmp.snduration(1:(iSn - 1)));
        end;
        onsets{n} = [onsets{n}; newonsets(:)];
        durations{n} = [durations{n}; newdurations(:)];
        if isfield(multi, 'pmod') && (numel(multi(1).pmod) >= n)
            for p = 1:numel(multi(1).pmod(n).param)
                pmod(n).param{p} = [pmod(n).param{p}; multi(iSn).pmod(n).param{p}(:)];
            end;
        end;
    end;
   
end;

% normalise if desired --
if cfg.norm
    Y = (Y - mean(Y))/std(Y);
end;
Y = Y(:);

% collect information into tmp --
tmp.length=numel(Y);

% scale pmods before orthogonalisation -- 
tmp.pmodno=zeros(numel(names), 1);
if exist('pmod', 'var')
    for n=1:numel(pmod)
        if ~isempty(pmod(n).param)
            for p=1:numel(pmod(n).param)
                % mean center and scale pmods
                try
                    tmp.pmodscale(n, p) = std(pmod(n).param{p});
                catch
                    tmp.pmodscale(n, p) = 1;
                end;
                pmod(n).param{p}=(pmod(n).param{p}-mean(pmod(n).param{p}))/tmp.pmodscale(n, p);
            end;
            % register number of pmods
            tmp.pmodno(n)=p;
        end;
    end;
else
    pmod = [];
end;

% collect data & regressors for output model --
glm.input.data    = y;
glm.input.sr      = sr;
glm.Y             = Y; 
glm.M             = M;
glm.infos.sr      = newsr;
glm.infos.duration     = numel(glm.Y)/glm.infos.sr;
glm.infos.durationinfo = 'duration in seconds';
glm.timing.multi      = multi;
glm.timing.names      = names;
glm.timing.onsets     = onsets;
glm.timing.durations  = durations;
glm.timing.pmod       = pmod;

% clear local variables --
clear iSn iMs ynew newonsets newdurations newmissing missingtimes 


% create temporary onset functions
%-------------------------------------------------------------------------
% cycle through conditions
for iCond = 1:numel(names)
    tmp.regscale{iCond} = 1;
    % first process event onset, then pmod
    tmp.onsets = onsets{iCond};
    tmp.durations = durations{iCond};
    % if file starts with first event, set that onset to 1 instead of 0
    if any(tmp.onsets == 0)
        tmp.onsets(tmp.onsets == 0) = 1;
    end;
    col=1;
    tmp.colnum=1+tmp.pmodno(iCond);
    tmp.X{iCond}=zeros(tmp.length, tmp.colnum);
    for k = 1:numel(tmp.onsets)
        tmp.X{iCond}(tmp.onsets(k):(tmp.onsets(k) + tmp.durations(k)), col)=1;
    end;
    tmp.name{iCond, col}=names{iCond};
    col=col+1;
    if exist('pmod') && ~isempty(pmod)
        if iCond<=numel(pmod)
            if ~isempty(pmod(iCond).param)
                for p=1:numel(pmod(iCond).param)
                    for k = 1:numel(tmp.onsets)
                        tmp.X{iCond}(tmp.onsets(k):(tmp.onsets(k) + tmp.durations(k)), col)=pmod(iCond).param{p}(k);
                    end;
                    tmp.name{iCond, col}=[names{iCond}, ' x ', pmod(iCond).name{p}];
                    tmp.regscale{iCond}(col) = tmp.pmodscale(iCond, col - 1);
                    col=col+1;
                end;
            end;
        end;
        % orthogonalize pmods before convolution
        foo = spm_orth(tmp.X{iCond});
        % catch zero matrices (unclear yet why this happens, 01-Apr-2012)
        if all(all(foo==0))
            warning('the pmods in condition %i have not been orthogonalized (because spm_orth returned a zero matrix)', iCond)
        else
            tmp.X{iCond} = foo;
        end
    end;
end;


% create design matrix
%-------------------------------------------------------------------------
% create design matrix filter
Xfilter = cfg.filter; 
Xfilter.sr = glm.infos.sr; 
Xfilter.down = 'none'; % turn off no low pass warning 

% convolve with basis functions
snoffsets = cumsum(tmp.snduration);
snonsets  = [1, snoffsets(1:end) + 1];
tmp.XC = cell(1,numel(names));
for iCond = 1:numel(names)
    tmp.XC{iCond} = [];
    tmp.regscalec{iCond} = [];
    iXCcol = 1;
    for iXcol = 1:size(tmp.X{iCond}, 2)
        for iBf = 1:glm.bf.bfno
%             figure;
%             plot(glm.bf.X(:,iBf),'Linewidth', 10)
%             grid on
%             set(gca,'fontsize',20)
%             
%             set(gca, 'xtick', 0:100:900)
%             oldxlabels=get(gca,'xtick')
%             newxlabels=oldxlabels/10
%             set(gca,'xticklabel',newxlabels)
%            
%          %  axis([0 900 0 1])
%           %   set(gca, 'xticklabel', ['10', '20', '30', '40', '50', '60', '70', '80', '90'])
% 
% 
%             xlabel('Time (sec)')
%             ylabel('Response function')
%             figure;
%             plot(glm.bf.X(1:150,iBf),'Linewidth', 10)
%             grid on
%             set(gca, 'ylim', [-0.2, 1.2])
            % process each session individually 
            for iSn = 1:numel(tmp.snduration)
                % convolve 
                tmp.col{iSn, 1} = conv(tmp.X{iCond}(snonsets(iSn):snoffsets(iSn), iXcol), glm.bf.X(:,iBf));
                % filter design matrix w/o downsampling
                [sts,  tmp.col{iSn, 1}] = scr_prepdata(tmp.col{iSn, 1}, Xfilter);
                if sts ~= 1, glm = struct([]); return; end;
                % cut away tail
                tmp.col{iSn, 1}((tmp.snduration(iSn) + 1):end) = [];
            end;
            tmp.XC{iCond}(:, iXCcol) = cell2mat(tmp.col);
            tmp.namec{iCond}{iXCcol, 1} = [tmp.name{iCond, iXcol}, ', bf ', num2str(iBf)];
            tmp.regscalec{iCond} = [tmp.regscalec{iCond}, tmp.regscale{iCond}(iXcol)];
            iXCcol = iXCcol + 1;
            % clear local variable
            tmp.col = {};
        end;
    end;
    
    % mean center
    for iXCol=1:size(tmp.XC{iCond},2)
        tmp.XC{iCond}(:,iXCol) = tmp.XC{iCond}(:,iXCol) - mean(tmp.XC{iCond}(:,iXCol));
    end;
    
    % orthogonalize after convolution if there is more than one column per
    % condition
    if size(tmp.XC{iCond}, 2) > 1
        foo=spm_orth(tmp.XC{iCond});
        % catch zero matrices (unclear yet why this happens, 01-Apr-2012)
        if ~(all(foo(:) == 0))
            tmp.XC{iCond} = foo;
        else
            warning('\nOrthogonalisation error in event type #%02.0f\nCorrelation coefficients are: ', iCond);
            cc = corrcoef(tmp.XC{iCond});
            for k = 2:size(cc, 1)
                fprintf('%0.2f, ', cc(1, k));
            end;
            fprintf('\n');
        end;
    end;
end;

% define model
glm.X = cell2mat(tmp.XC);
glm.regscale = cell2mat(tmp.regscalec);
r=1;
for iCond = 1:numel(names)
    n = numel(tmp.namec{iCond});
    glm.names(r:(r+n-1), 1) = tmp.namec{iCond};
    r = r + n;
end;

% add nuisance regressors
for iSn = 1:numel(cfg.datafile)
    Rf{iSn} = [];
    cfg.filter.sr = sr(iSn);
    for iR = 1:nR
        [sts, Rf{iSn}(:, iR)]  = scr_prepdata(R{iSn}(:, iR), cfg.filter);
        if sts ~= 1, return; end;
    end
end
Rf = cell2mat(Rf(:));

for iR = 1:nR
    glm.names{end+1, 1} = sprintf('R%01.0f', iR);
end;

glm.X = [glm.X, Rf];
glm.regscale((end+1):(end+nR)) = 1;

% add constant(s)
r=1;
for iSn = 1:numel(cfg.datafile);
    glm.X(r:(r+tmp.snduration(iSn)-1), end+1)=1;
    glm.names{end+1, 1} = ['Constant ', num2str(iSn)];
    r = r + tmp.snduration(iSn);
end;
glm.interceptno = iSn;
glm.regscale((end+1):(end+iSn)) = 1;

% delete missing epochs and prepare output
glm.YM = glm.Y;
glm.YM(glm.M==1) = [];
glm.Y(glm.M==1) = NaN;
glm.XM = glm.X;
glm.XM(glm.M==1, :) = [];
glm.X(glm.M==1, :) = NaN;
glm.Yhat    = NaN(size(Y));

% clear local variables
clear tmp Xfilter r iSn n iCond


% invert model & save
%-------------------------------------------------------------------------
% this is where the beef is
glm.stats = pinv(glm.XM)*glm.YM;           % parameter estimates
glm.Yhat(glm.M==0) = glm.XM*glm.stats;     % predicted response
glm.e    = glm.Y - glm.Yhat;               % residual error
glm.EV   = 1 - (var(glm.e)/var(glm.YM));   % explained variance proportion

% rescale pmod parameter estimates & design matrix
%-------------------------------------------------------------------------
glm.X = glm.X .* repmat(glm.regscale, size(glm.X, 1), 1);
glm.XM = glm.XM .* repmat(glm.regscale, size(glm.XM, 1), 1);
glm.stats = glm.stats .* glm.regscale';

savedata = struct('glm', glm);
savedata.modeltype = 'glm';
scr_load1(cfg.modelfile, 'save', savedata, options);



[sts, glm] = scr_load1(cfg.modelfile, 'all', 'glm');
if sts ~= 1, return; end;
bs = glm.bf.X;
bfno = glm.bf.bfno;
bfdur = size(bs, 1);

% find all non-nuisance regressors 
% -------------------------------------------------------------------------
regno = (numel(glm.names) - glm.interceptno)/bfno;
if regno ~= floor(regno), warning('Mismatch in basis functions and regressor number.'); return; end;

% reconstruct responses
% -------------------------------------------------------------------------
resp = NaN(bfdur, regno);
for k = 1:regno
    condname{k} = glm.names{((k - 1) * bfno + 1)};
    foo = strfind(condname{k}, ', bf');
    condname{k} = [condname{k}(1:(foo-1)), ' recon']; clear foo;
    resp(:, k) = bs * glm.stats(((k - 1) * bfno + 1):(k * bfno));
    [foo, bar] = max(abs(resp(:, k)));
    recon(k, 1) = resp(bar, k);
end;

% save
% -------------------------------------------------------------------------
glm.recon = recon;
glm.resp  = resp;
glm.reconnames = condname(:);
save(cfg.modelfile, 'glm');
modelfile=cfg.modelfile;

% user output
%-------------------------------------------------------------------------

fprintf(' done. \n');

return

