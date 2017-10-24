%%
function [ score ] = test_fir_pop_member( currMem, targInput, makePlot )

if nargin<3, makePlot = false; end

%% Initialize values that have been commanded
rip = targInput.rip; % Ripple in dB
rej = targInput.rej; % Rejection in dB
fpass = targInput.fpass; % PB corner frequency
fstop = targInput.fstop; % SB frequency
pbscale = targInput.pbscale; % Scale of PB score
sbscale = targInput.sbscale; % Scale of SB score
fs = targInput.fs; % sample freq
nTap = targInput.nTap; % number of coeffs
fftlen = targInput.fftlen; % fftlength to test
sbramp = targInput.sbramp; % stop-band ramp function
ScrFnc = targInput.ScrFnc;
if(isfield(targInput,'isScl')), isScl = targInput.isScl; else, isScl = 1;end
if(isfield(targInput,'osScl')), osScl = targInput.osScl; else, osScl = 1;end

%% Generate the FFT of our member
% Test whether we're making an even or odd length filter
if 2*round(nTap/2)==nTap
    b_a = [currMem, fliplr(currMem)] ;
else
    b_a = [currMem, fliplr(currMem(1:end-1))] ;
end
% Normalize the member under test before comparing
b_a = b_a./sum(b_a);
% Take the FFT to compare to input filter
[h_a, f] = freqz(b_a, 1, fftlen, fs);
ha_dB = 20*log10(abs(h_a));

%%
% This is just the filter we want to match for reference - it has nothing
% to do with our scoring
h_f = targInput.h_targ;
hf_dB = 20*log10(abs(h_f));

%% Score the PB of the filter
passix = (f<=fpass); % Boolean vector correlating to the bins in our passband
% PB comparison
compPb = hf_dB(passix);
meancpb = median(compPb);
% Trim off the PB data for scoring
apb = ha_dB(passix); % We compare PB as a dB value
% Find the median of the PB data - this is what we'll be scoring PB against
meanpb = median(apb);
% Find the error in our PB
pberr = abs(apb - meanpb)./rip;
% Generate a scaling vector based on whether we're in or out of ripple
pbScl = pberr<=rip; pbScl(pbScl==1) = isScl; pbScl(pbScl==0) = osScl;

%% Score the SB of the filter
stopix = (f>=fstop); % Boolean vector correlating to the bins in our stopband
rmpScl = linspace(sbramp(1),sbramp(2),sum(stopix)).';
% SB Comparison
compSb = hf_dB(stopix); % This is to check against our rejection requirement
% Trim off the SB data for scoring
asb = ha_dB(stopix); % This is to check against our rejection requirement
linsb = h_a(stopix); % We want to look at linear SB to get as close to 0 as possible
% Find the SB value to compare to, generate a scaling vector based on that
sbSpec = meanpb - rej;
sbLinSpec = 10^(sbSpec/20);
sbScl = asb<=sbSpec; sbScl(sbScl==1) = isScl; sbScl(sbScl==0) = osScl;
% The stop band error is just equal to the stob-band value, we want to get
% as close to 0 as possible
sberr = abs(linsb)./sbLinSpec;

%% Max Errors
maxPbErr = 0;%max(pberr);
maxSbErr = 0;%max(sberr);

%% Compile the overall score
% Now score the PB and SB
% pbscore = pbscale*(mean((pberr.*pbScl).^4) + maxPbErr^4);
% sbscore = sbscale*(mean(rmpScl.*(sberr.*sbScl).^4) + maxSbErr^4);
pbscore = pbscale*( ScrFnc(pberr.*pbScl) + ScrFnc(maxPbErr) );
sbscore = sbscale*( ScrFnc(rmpScl.*sberr.*sbScl) + ScrFnc(maxSbErr) );
score = pbscore + sbscore;

%% Plot the results
if makePlot
    fminp = min(f(passix));
    fmaxp = max(f(passix));
    
    figure(202);
    
    subplot(221);
    plot( f(passix), compPb - meancpb, '-m', 'LineWidth', 2.0 );hold on;
    plot( f(passix), apb - meanpb, 'LineWidth', 2.0 );
    plot([fminp, fmaxp], [0,0], '--m', 'LineWidth', 2.0 );
    plot([fminp, fmaxp], [rip,rip], '--m', 'LineWidth', 2.0 );
    plot([fminp, fmaxp], -[rip,rip], '--m', 'LineWidth', 2.0 );hold off
    grid;ylabel('Mag.[dB]'); zoom on;
    title(sprintf('Passband Response (Fp = %0.2f MHz)\nNtaps=%g, pbscore=%g', fpass, nTap, pbscore ) );
    if(rip==0),ylim([-0.1,0.1]),else, ylim(2.*[-rip,rip]), end
    
    subplot(222);
    plot( f(stopix), compSb - meancpb, '-m', 'LineWidth', 2.0 );hold on;
    plot( f(stopix), asb - meanpb, 'LineWidth', 2.0 );hold on;
    plot([min(f(stopix)), max(f(stopix))], [-rej,-rej], '--m', 'LineWidth', 2.0 );hold off
    grid; zoom on;
    title(sprintf('Stopband Response (Fs = %0.2f MHz)\nNtaps=%g, sbscore=%g', fstop, nTap, sbscore ) );
    ylim(-rej+[-10,30])
    
    subplot(212);
    plot( f, hf_dB - meancpb, '-m','LineWidth', 2.0 ); hold on
    plot( f, ha_dB - meanpb,'LineWidth', 2.0 ); hold off
    grid; xlabel('Freq[MHz]');ylabel('Mag.[dB]');zoom on;
    title(sprintf('Magnitude Response, Ntaps=%g, score=%g', nTap, score ) );
    ylim([-100,3])
    drawnow;
    
end





end
