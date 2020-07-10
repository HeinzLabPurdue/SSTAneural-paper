% function [sig, time]= create_SAM(fc, fm, fs, modDepth, dur, amp, phi_c, phi_m)
function [sig, time]= create_SAM(fc, fm, fs, modDepth, dur, amp, phi_c, phi_m)

if nargin <3
    error ('Need more love');
end

if ~exist('phi_m', 'var')
    phi_m= -pi/2;
elseif isempty(phi_m)
    phi_m= -pi/2;
end

if ~exist('phi_c', 'var')
    phi_c= 0;
elseif isempty(phi_c)
    phi_c= 0;
end

if ~exist('amp', 'var')
    amp= 1;
elseif isempty(amp)
    amp= 1;
end

if ~exist('dur', 'var')
    dur= .2;
elseif isempty(dur)
    dur= .2;
end

if ~exist('modDepth', 'var')
    modDepth= 1;
elseif isempty(modDepth)
    modDepth= 1;
end


time= (0:1/fs:dur-1/fs)';
sig= amp*(1+modDepth*sin(2*pi*fm*time + phi_m)).*sin(2*pi*fc*time + phi_c);
