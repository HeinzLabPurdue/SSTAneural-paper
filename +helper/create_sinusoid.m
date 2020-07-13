% function [sig, time]=create_sinusoid(freq, fs, dur, amp, phi)
function [sig, time]=create_sinusoid(freq, fs, dur, amp, phi)

if nargin<2
    error('need atleast two input');
end

if ~exist('dur', 'var')
    dur= 1;
elseif isempty(dur)
    dur= 1;
end

if ~exist('amp', 'var')
    amp= ones(size(freq));
elseif isempty(amp)
    amp= ones(size(freq));
elseif (numel(freq)~= numel(amp)) && (numel(amp) == 1)
    amp= amp*ones(size(freq));
elseif numel(amp) ~= numel(freq)
    error('amplitude and frequency are of different length');
end

if ~exist('phi', 'var')
    phi= 2*pi*randn(size(freq));
elseif isempty(phi)
    phi= 2*pi*randn(size(freq));
elseif (numel(freq)~= numel(phi)) && (numel(phi) == 1)
    phi= phi*ones(size(freq));
elseif numel(phi) ~= numel(freq)
    error('phase and frequency are of different length');
end


time= (0:1/fs:dur-1/fs)';
sig=zeros(size(time));

for freqVar=1:length(freq)
   sig= sig+ amp(freqVar)*sin(2*pi*freq(freqVar)*time + phi(freqVar));
end