% Copyright (C) 2017 Trevor
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @deftypefn {} {@var{retval} =} TestPopulation (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: Trevor <Trevor@TREVOR>
% Created: 2017-10-21

function [score, bstIx] = TestPopulation(fTaps, targInput, makePlot)

%% Initialize inputs
if(nargin<3), makePlot = false; end

%% Initialize values that have been commanded
rip       = targInput.rip;     % Ripple in dB
rej       = targInput.rej;     % Rejection in dB
fpass     = targInput.fpass;   % PB corner frequency
fstop     = targInput.fstop;   % SB frequency
pbscale   = targInput.pbscale; % Scale of PB score
tbscale   = targInput.tbscale;
sbscale   = targInput.sbscale; % Scale of SB score
fs        = targInput.fs;      % Sample freq
fftlen    = targInput.fftlen;  % Fftlength to test
ScrFnc    = targInput.ScrFnc;  % Function used to score filter

%% Generate the FFT for comparison
f_full    = (0:1/fftlen:1-1/fftlen).*fs;           % Generate a frequency vector
f         = f_full(1:fftlen/2);                    % Trim off the unique values
h_full    = fft(fTaps.'./sum(fTaps.'), fftlen).';  % Generate the FFT for comparison
h         = h_full(:, 1:fftlen/2);                 % Trim off the unique values
h_dB      = 20*log10(abs(h));                      % Convert to dB for scoring

%% Score the PB of the filter
passix    = (f <= fpass);             % Boolean vector correlating to the bins in our passband
apb       = h_dB(:, passix);          % Trim off the PB data for scoring, we compare PB as a dB value
medianpb  = median(apb, 2);           % Find the median of the PB data - this is what we'll be scoring PB against
pberr     = abs(apb - medianpb)./rip; % Find the error in our PB relative to the ripple

%% Score the SB of the filter
stopix    = (f >= fstop);             % Boolean vector correlating to the bins in our stopband
linsb     = h(:, stopix);             % We want to look at linear SB to get as close to 0 as possible
sbSpec    = medianpb - rej;           % Find the SB value to compare to, generate a scaling vector based on that
sbLinSpec = 10.^(sbSpec/20);          % Get the linear spec for scoring reference
sberr     = abs(linsb)./sbLinSpec;    % Get the SB as close to 0 as possible

%% Transition band - make sure it decreases
atb = h(:, (f>fpass)&(f<fstop) );
tberr = abs(atb(:, 2:end))./abs(atb(:, 1:end-1));
% tberr = (abs(atb(:, 2:end))./abs(atb(:, 1:end-1)))<1;

%% Compile the overall score
pbscore = pbscale*ScrFnc(pberr);      % Score the PB
tbscore = tbscale*ScrFnc(tberr);      % Score the SB
sbscore = sbscale*ScrFnc(sberr);      % Score the SB
score   = pbscore + tbscore + sbscore;          % Total score
[bstScore, bstIx] = min(score);

%% Plot the results
if makePlot
    fminp             = min(f(passix));
    fmaxp             = max(f(passix));
    asb               = h_dB(:, stopix);     % Trim off the SB data for scoring, this is to check against our rejection requirement
    bstPbScore        = pbscore(bstIx);
    bstTbScore        = tbscore(bstIx);
    bstSbScore        = sbscore(bstIx);
    figure(202);
    
    subplot(221);
    plot( f(passix), apb(bstIx,:) - medianpb(bstIx,:), 'LineWidth', 2.0 );hold on;
    plot([fminp, fmaxp], [0,0], '--m', 'LineWidth', 2.0 );
    plot([fminp, fmaxp], [rip,rip], '--m', 'LineWidth', 2.0 );
    plot([fminp, fmaxp], -[rip,rip], '--m', 'LineWidth', 2.0 );hold off
    grid;ylabel('Mag.[dB]'); zoom on;
    title(sprintf('Passband Response (%0.2f MHz)\nBest PB Score=%g', fpass, bstPbScore ) );
    if(rip==0),ylim([-0.1,0.1]),else, ylim(2.*[-rip,rip]), end
    
    subplot(222);
    plot( f(stopix), asb(bstIx,:) - medianpb(bstIx,:), 'LineWidth', 2.0 );hold on;hold on;
    plot([min(f(stopix)), max(f(stopix))], [-rej,-rej], '--m', 'LineWidth', 2.0 );hold off
    grid; zoom on;
    title(sprintf('Stopband Response (%0.2f MHz)\nBest SB Score=%g', fstop, bstSbScore ) );
    ylim(-rej+[-10,30])
    
    subplot(212);
    plot( f, h_dB(bstIx,:) - medianpb(bstIx,:),'LineWidth', 2.0 ); hold on
    plot([fminp, fmaxp], [rip,rip], '--m', 'LineWidth', 2.0 );
    plot([fminp, fmaxp], -[rip,rip], '--m', 'LineWidth', 2.0 );
    plot([min(f(stopix)), max(f(stopix))], [-rej,-rej], '--m', 'LineWidth', 2.0 );hold off
    grid; xlabel('Freq[MHz]');ylabel('Mag.[dB]');zoom on;
    title(sprintf('Magnitude Response Best, TB Score=%g\nBest Score=%g', bstTbScore, bstScore ) );
    ylim([-100,3])
    set(figure(202),'Position',[300 200 560 400 ]);
    drawnow;
    
end




end
