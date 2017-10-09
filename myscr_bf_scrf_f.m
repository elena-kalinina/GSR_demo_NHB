
function [bs] = myscr_bf_scrf(varargin)

%EXAMPLE CALL
%mybs=myscr_bf_scrf(1, 2);

% SCR_infbs constructs an informed basis set with a biexponentially 
% modified gaussian function and derivatives to time and dispersion
%
% FORMAT: [INFBS] = SCR_BF_SCRF(TD, D) 
%     OR: [INFBS] = SCR_BF_SCRF([TD, D]) 
% with td = time resolution in s and d:number of derivatives (default 0)
%
% REFERENCE
% Bach DR, Flandin G, Friston KJ, Dolan RJ (2010). Modelling event-related 
% skin conductance responses. International Journal of Psychophysiology,
% 75, 349-356.
%________________________________________________________________________
% PsPM 3.0
% (C) 2008-2015 Dominik R Bach (Wellcome Trust Centre for Neuroimaging)

% $Id: scr_bf_scrf.m 702 2015-01-22 15:06:14Z tmoser $   
% $Rev: 702 $

% initialise
global settings
if isempty(settings), scr_init; end;

% check input arguments
if nargin==0
    errmsg='No sampling interval stated'; warning(errmsg); return;
end;

td = varargin{1}(1);
if numel(varargin{1}) == 1 && nargin == 1
    d = 0;
elseif numel(varargin{1}) == 2
    d = varargin{1}(2);
else
    d = varargin{2}(1);
end;
    
if (d<0)||(d>2), d=0; end;

% get parameters and basis function
[bs(:, 1), p] = myscr_bf_scrf_f(td);
if d>0
    bs(:, 2) = [0; diff(bs(:,1))]; 
    bs(:, 2) = bs(:,2)/sum(abs(bs(:,2)));
    figure;
    plot(bs(:, 2));
end;
if d>1 
    p(2) = 1.8 * p(2); 
    bs(:, 3) = bs(:, 1) - myscr_bf_scrf_f(td, p); 
    bs(:, 3) = bs(:, 3)/sum(abs(bs(:, 3)));
    figure;
    plot(bs(:, 3));
    
end;

% orthogonalize
bs=spm_orth(bs);
return;


function [ft, p] = myscr_bf_scrf_f(td, p)
% SCR_bf_crf: canonical skin conductance response function 
% (exponentially modified gaussian, EMG)
% FORMAT: [bf p] = SCR_bf_scrf_f(td, p)
% with  td = time resolution in s
%       p(1): time to peak
%       p(2): variance of rise defining gaussian
%       P(3:4): decay constants
% 
% REFERENCE
% Bach DR, Flandin G, Friston KJ, Dolan RJ (2010). Modelling event-related skin 
% conductance responses. International Journal of Psychophysiology, 75, 349-356.
%__________________________________________________________________________
% PsPM 3.0
% (C) 2009-2015 Dominik R Bach (Wellcome Trust Centre for Neuroimaging)

% $Id: scr_bf_scrf_f.m 702 2015-01-22 15:06:14Z tmoser $   
% $Rev: 702 $

% initialise
% -------------------------------------------------------------------------
global settings;
if isempty(settings), scr_init; end;
% -------------------------------------------------------------------------

if nargin < 1
   errmsg='No sampling interval stated'; warning(errmsg); return;
elseif nargin < 2
    p=[3.0745  0.7013 0.3176 0.0708];
end;

t0 = p(1);
sigma = p(2);
lambda1 = p(3);
lambda2 = p(4);

t = (td:td:90)';

gt = exp(-((t - t0).^2)./(2.*sigma.^2));
ht = exp(-t*lambda1) + exp(-t*lambda2);

ft = conv(gt, ht);
ft = ft(1:numel(t));
ft = ft/max(ft);

% figure;
% plot(gt);
% figure;
% plot(ht);
% figure;
% plot(ft);
