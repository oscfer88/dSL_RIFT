%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function [XDL, XDR, pDelays] = cross_correlation(X, DL, DR, tTrialBeg, tTrialEnd)

% trial info
%tTrialBeg = -0.2;
%tTrialEnd = 1.7;
bCAT = 1; % concatenate over trils (1, default) or not (0)

nSamples = size(X, 2);
nTrials = size(X, 3);
t = linspace(tTrialBeg, tTrialEnd, nSamples);
iTrialBeg = find(t >= 0, 1);

% MEG % tagging parameters
fs = 1000; % sampling rate
lSegment = 1.5 * fs; % segment for correlation (in seconds)
fl = 45; % low tagging frequency
fh = 85; % high tagging frequency

% set delay range
pDelays = -200:200;
nDelays = length(pDelays);

DL = squeeze(DL);
DR = squeeze(DR);

% filter data 
for iTrial = 1:nTrials
  y = X(:, :, iTrial); y = apply_filter(y, fl, fh, fs); X(:, :, iTrial) = y; 
end

% filter tagging signal 
y = apply_filter(DL', fl, fh, fs); DL = y'; 
y = apply_filter(DR', fl, fh, fs); DR = y';

% flatten data
if bCAT == 1 
  X = X(:, iTrialBeg:(iTrialBeg + lSegment - 1), :);
  DL = DL(iTrialBeg:(iTrialBeg + lSegment - 1), :);
  DR = DR(iTrialBeg:(iTrialBeg + lSegment - 1), :);
  X = reshape(X, size(X, 1), []);
  DL = DL(:);
  DR = DR(:);
  lSegment = length(DL) - 2 * length(pDelays);
end

% compute xcorr
nChannels = size(X, 1);
nTrials = size(X, 3);

% init
XDL = zeros(nChannels, nDelays, nTrials); 
XDR = zeros(nChannels, nDelays, nTrials); 

% over trials
iBeg = find(t >= 0, 1);
% if bCAT == 1
%     iBeg = 1;
% end
iEnd = iBeg + lSegment - 1;
for iTrial = 1:nTrials
  yDL = DL(iBeg:iEnd, iTrial);
  yDR = DR(iBeg:iEnd, iTrial);
  yDL = repmat(yDL', nChannels, 1);
  yDR = repmat(yDR', nChannels, 1);
  % over delays
  for iDelay = 1:nDelays
    fprintf(1, '%d/%d | %d/%d\n', iTrial, nTrials, iDelay, nDelays);
    d = pDelays(iDelay);
    yX = X(:, (iBeg + d):(iEnd + d), iTrial);
    XDL(:, iDelay, iTrial) = diag(corr(yX', yDL'));
    XDR(:, iDelay, iTrial) = diag(corr(yX', yDR')); 
  end
end

end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function X = apply_filter(X, fl, fh, fs)

nPadLen = floor(0.2 * size(X, 2)); % 20%

% high pass filter
[b, a] = butter(4, fl / (fs / 2), 'high'); 
X = filtfilt(b, a, [X(:, nPadLen:-1:1), X, X(:, end:-1:(end - nPadLen))]')'; 
X(:, (end - nPadLen):end) = [];
X(:, 1:nPadLen) = [];

% low pass filter
[b, a] = butter(4, fh / (fs / 2), 'low'); 
X = filtfilt(b, a, [X(:, nPadLen:-1:1), X, X(:, end:-1:(end - nPadLen))]')'; 
X(:, (end - nPadLen):end) = [];
X(:, 1:nPadLen) = [];

end % end

%-------------------------------------------------------------------------------

