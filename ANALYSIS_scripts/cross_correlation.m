%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function cross_correlation(X, y, nTrials)

% MEG % tagging parameters
fs = 1000; % sampling rate
fl = 55; % low tagging freqeuncy
fh = 75; % high tagging freqeuncy

% filter data and tagging signal
for iTrial = 1:nTrials
  X(:, :, iTrial) = apply_filter(X(:, :, iTrial), fl, fh, fs);
  y(:, :, iTrial) = apply_filter(y(:, :, iTrial), fl, fh, fs);
  y(:, :, iTrial) = y(:, :, iTrial) ./ abs(hilbert(y(:, :, iTrial)')');
end

% concatenate over trials
X = reshape(X, size(X, 1), []);
y = reshape(y, size(y, 1), []);

% cross-correlation
pDelays = -200:1:200;
A = get_peak_delay(X', y', pDelays);

% plot
plot(pDelays, A); xlabel('delays'); ylabel('correlation'); legend('MEG1922','MEG1923','MEG2032','MEG2033','MEG2042','MEG2043','MEG2112','MEG2113','MEG2122','MEG2123','MEG2342','MEG2343');
figure; plot(pDelays, abs(hilbert(A'))); xlabel('delays'); ylabel('correlation'); legend('MEG1923','MEG1922','MEG2032','MEG2033','MEG2042','MEG2043','MEG2113','MEG2112','MEG2122','MEG2123','MEG2343','MEG2342');

end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function A = get_peak_delay(X, y, pDelays)

% autocorrelation
nChannels = size(X, 2);
L = length(pDelays);
A = zeros(nChannels, L);
for i = 1:length(pDelays)
  p = pDelays(i);
  if p == 0
    p = 1;
  end
  if p > 0
    u = [y(p:end); y(1:(p - 1))];
  else
    u = [y((end + p + 1):end); y(1:(end + p))];
  end
  U = repmat(u, 1, nChannels);
  % COR
  t = diag(corr(X, U));
  A(:, i) = t;
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

