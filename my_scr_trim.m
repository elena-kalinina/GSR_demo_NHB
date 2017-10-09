function newdatafile = my_scr_trim(datafile, cfg)
% SCR_TRIM cuts an SCR dataset to the limits set with the parameters 'from'
% and 'to' and writes it to a file with a prepended 't'
%
% FORMAT:
% NEWDATAFILE = SCR_TRIM (DATAFILE, FROM, TO, REFERENCE, options)
%
% datafile: a file name, a cell array of filenames, a struct with
%           fields .data and .infos or a cell array of structs
% from and to: either numbers, or 'none'
% reference: 'marker': from and to are set in seconds with 
%                         respect to the first and last scanner/marker pulse
%            'file':    from and to are set in seconds with respect to start
%                         of datafile 
%            a 2-element vector: from and to are set in seconds with
%                         respect to the two markers defined here
%
% options:  options.overwrite:       overwrite existing files by default
%           options.marker_chan_num: marker channel number - if undefined 
%                                     or 0, first marker channel is used
%
% RETURNS a filename for the updated file, a cell array of filenames, a
% struct with fields .data and .infos or a cell array of structs
%
%__________________________________________________________________________
% PsPM 3.0
% (C) 2008-2015 Dominik R Bach (Wellcome Trust Centre for Neuroimaging)

% $Id: scr_trim.m 701 2015-01-22 14:36:13Z tmoser $
% $Rev: 701 $

% initialise
% -------------------------------------------------------------------------
global settings;
if isempty(settings), scr_init; end;
newdatafile = [];

% check input arguments
% -------------------------------------------------------------------------
% if nargin<1
%     warning('ID:invalid_input', 'No data.\n'); return;
% elseif nargin<2
%     warning('ID:invalid_input', 'No start or end point given.\n'); return;
% elseif nargin<3
%     warning('ID:invalid_input', 'No end point given.\n'); return;
% elseif nargin<4
%     warning('ID:invalid_input', 'No reference given.\n'); return;
% end;
% 
% if ~((ischar(from) && strcmpi(from, 'none')) || (isnumeric(from) && numel(from) == 1)) 
%     warning('ID:invalid_input', 'No valid start point given.\n'); return;
% elseif ~((ischar(to) && strcmpi(to, 'none')) || (isnumeric(to) && numel(to) == 1))
%     warning('ID:invalid_input', 'No end point given'); return;
% end;

if nargin <2
    cfg=[];
end

if ~isfield(cfg, 'reference')
    
    cfg.reference='file';
       
end

if ~isfield(cfg, 'from')
    
    cfg.from=120;
       
end

if ~isfield(cfg, 'to')
    
    cfg.to=360;
       
end

if ~isfield(cfg, 'options')
    
    cfg.options=[];
    cfg.options.overwrite=0; %1 for overwrite
    
end




if strcmpi(cfg.reference, 'marker')
    getmarker = 1;
    startmarker = 1;
    endmarker = [];
elseif isnumeric(cfg.reference) && numel(cfg.reference) == 2
    getmarker = 1;
    startmarker = cfg.reference(1);
    endmarker = cfg.reference(2);    
    cfg.reference = 'events';
    % check if reference markers are valid ---
    if startmarker < 1 || endmarker < startmarker
        warning('ID:invalid_input', 'No valid reference markers.\n'); return;
    end;
elseif strcmpi(cfg.reference, 'file')
    getmarker = 0;
else
    warning('ID:invalid_input', 'Invalid reference option ''%s''', cfg.reference); return;
end;

% set options ---
try
    cfg.options.overwrite; 
catch
    cfg.options.overwrite = 0;
end;

if ~isfield(cfg.options,'marker_chan_num') || ~isnumeric(cfg.options.marker_chan_num) || numel(cfg.options.marker_chan_num) > 1
    cfg.options.marker_chan_num = 0;
end;

% check data file argument --
if ischar(datafile) || isstruct(datafile)
    D = {datafile};
elseif iscell(datafile) 
    D = datafile;
else
    warning('Data file must be a char, cell, or struct.');
end;
clear datafile

% work on all data files
% -------------------------------------------------------------------------
for d=1:numel(D)
    % determine file names ---
    datafile=D{d};
        
    % user output ---
    if isstruct(datafile)
        fprintf('Trimming ... ');
    else
        fprintf('Trimming %s ... ', datafile);
    end;
    
    % check and get datafile ---
    [sts, infos, data] = scr_load_data(datafile, 0);
    if getmarker
        if cfg.options.marker_chan_num
            [nsts, ninfos, ndata] = scr_load_data(datafile, cfg.options.marker_chan_num);
            if ~strcmp(ndata{1}.header.chantype, 'marker')
                warning('ID:invalid_option', 'Channel %i is no marker channel. The first marker channel in the file is used instead', cfg.options.marker_chan_num);
                [nsts, ninfos, ndata] = scr_load_data(datafile, 'marker');
            end
        else
            [nsts, ninfos, ndata] = scr_load_data(datafile, 'marker');
        end;
        sts = [sts; nsts];
        events = ndata{1}.data;
        if isempty(endmarker), endmarker = numel(events); end;
        clear nsts ninfos ndata
    end;
    if any(sts == -1), newdatafile = []; break; end;
    
    % convert from and to into time in seconds ---
    if ischar(cfg.from) % 'none'
        startpoint = 0;
    else
        if getmarker % 'marker'
            startpoint = events(startmarker) + cfg.from;
        else         % 'file'
            startpoint = cfg.from;
        end;
    end;
    if ischar(cfg.to) % 'none'
        endpoint=infos.duration;
    else
        if getmarker  % 'marker'
            if endmarker > numel(events)
                warning('ID:marker_out_of_range', '\nEnd marker (%03.0f) out of file - no trimming end end.\n', endmarker);
                endpoint = infos.duration;
            else
                endpoint=events(endmarker) + cfg.to;
            end;
        else          % 'file'
            endpoint = cfg.to;
        end;
    end;
    
    % check start and end points ---
    if (startpoint < 0)
        warning('ID:marker_out_of_range', '\nStart point (%.2f s) outside file, no trimming at start.', startpoint);
        startpoint = 0;
    end;
    if endpoint > infos.duration
        warning('ID:marker_out_of_range', '\nEnd point (%.2f s) outside file, no trimming at end.', endpoint);
        endpoint = round(length(data{1}.data)/cfg.rec_sr);    %infos.duration;
    end;
    
    % trim file ---
    for k = 1:numel(data)
        if ~strcmpi(data{k}.header.units, 'events') % waveform channels
            % set start and endpoint
            newstartpoint = floor(startpoint * cfg.rec_sr); %data{k}.header.sr);
            if newstartpoint == 0, newstartpoint = 1; end;
            newendpoint = floor(endpoint * cfg.rec_sr); %data{k}.header.sr);
            if newendpoint > numel(data{k}.data), newendpoint = numel(data{k}.data); end;
            % trim data
            data{k}.data=data{k}.data(newstartpoint:newendpoint);
        else                                        % event channels
            data{k}.data(data{k}.data > endpoint) = [];
            data{k}.data = data{k}.data - startpoint;
            data{k}.data(data{k}.data < 0) = [];
        end;
        % save new file
        infos.duration = endpoint - startpoint;
        infos.trimdate = date;
        infos.trimpoints = [startpoint endpoint];
    end;
    clear savedata
    savedata.data = data; savedata.infos = infos; 
    if isstruct(datafile)
        sts = scr_load_data(savedata, 'none');
        newdatafile = savedata;
    else
        [pth, fn, ext] = fileparts(datafile);
        newdatafile    = fullfile(pth, ['t', fn, ext]);
        savedata.infos.trimfile = newdatafile;
        savedata.options = cfg.options;
        sts = my_scr_load_data(newdatafile, savedata);
    end;
    if sts ~= 1
        warning('Trimming unsuccessful for file %s.\n', newdatafile); 
    else
        Dout{d} = newdatafile;
        % user output
        fprintf('  done.\n');
    end;
end;

% if cell array of datafiles is being processed, return cell array of
% filenames
if d > 1
    clear newdatafile
    newdatafile = Dout;
end;

return;
