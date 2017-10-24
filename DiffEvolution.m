%% Differential Evolution Function

% Control variables for the cost function
close all force
clear variables;
format long g;

%% Top Level Controls
normPop  = true;                               % Normalize and make fixed-point the members
npf      = 10;                                 % # of pops
nTap     = 64;                                 % Number of filter coefficients
oddOrder = 2*round(nTap/2)==nTap;              % Determine if we have an odd order filter
refVal   = 50;
refScore = 200;

%% Differential Evolution Controls
strategy = 4;
itmax    = 1e6;                                % Maximum number of iterations to run
minScore = 1e-6;                               % Minimum Score to Achieve
D        = nTap/2;                             % Parameter Length
NP       = npf*D;                              % Population Size
F1       = 0.8;                                % Mutation Constant [0, 2]
F2       = 0.8;                                % Mutation Constant [0, 2]
CR       = 0.75;                               % Crossover Rate [0, 1]?
vMax     = 1;                                  % Upper Bound
vMin     = -1;                                 % Lower Bound

%% Filter Optimization Configurations
n                 = 2;                         % Exponent to use in cost function
targInput.ScrFnc  = @(x) (mean(abs(x).^n, 2)); % Cost function to use
targInput.rip     = 0.006;                     % pband peak variation in dB
targInput.rej     = 70;                        % dB rejection
targInput.fpass   = 8.1;                       % Passband corner freq in MHz
targInput.fstop   = targInput.fpass*1.85;      % targInput.fpass*9.175/4.865; % Stopband freq in MHz
targInput.pbscale = 10;                        % Scale factor for PB
targInput.tbscale = 10;                        % Scale factor for SB
targInput.sbscale = 1;                         % Scale factor for SB
targInput.fs      = 100;                       % Reference sample frequency MHz
targInput.nTap    = 64;                        % Number of filter coeffs
targInput.fftlen  = 2^10;                      % fft length for scoring
targInput.fpass   =  ceil( targInput.fpass/... % Quantization of fpass and fstop to the fft resolution
    (targInput.fs/targInput.fftlen/2) )*(targInput.fs/targInput.fftlen/2);
targInput.fstop   = floor( targInput.fstop/... % Quantization of fpass and fstop to the fft resolution
    (targInput.fs/targInput.fftlen/2) )*(targInput.fs/targInput.fftlen/2);

%% Start Differential Evolution Function
pop   = vMin + (vMax-vMin)*rand(NP, D);        % Initialize the population
pop = pop./(2*sum(pop,2));
pop = round( pop.*2.^( 15 - ceil(log2(max( abs(pop), [], 2 ))) ) );

%% Initialization for the DE algorithm
[score, bstIx] = TestPopulation([pop, fliplr(pop)], targInput, true); % Test the initial population
bm = pop(bstIx,:);
a = zeros(NP, 5);                              % Pre allocate random index matrix
cnt = 1;
compVec = 1:D;
rot = (0:1:NP-1);               % rotating index array (size NP)

while cnt<itmax || score<minScore
    popold = pop;
    
    for ix = 1:NP
        ind          = randperm(NP);   % Generate a random vector of choices
        ind(ind==ix) = [];             % Remove the current index from choice
        a(ix,:)      = ind(1:5);       % Take the first 5 random values
    end
    
    x1 = pop(a(:,1),:);                           % shuffled population 1
    x2 = pop(a(:,2),:);                           % shuffled population 2
    x3 = pop(a(:,3),:);                           % shuffled population 3
    x4 = pop(a(:,4),:);                           % shuffled population 4
    x5 = pop(a(:,5),:);                           % shuffled population 5
    
    %% Mutation Using DE Strategy
    switch strategy
        case 1 % DE/rand/1
            vi = x1 + F1*(x2 - x3);
        case 2 % DE/best/1
            vi = bm + F1*(x2 - x3);
        case 3 % DE/rand to best/1
            vi = x1 + F1*(x2 - x3) + F2*(bm - x1);
        case 4 % DE/curr. to best/1
            vi = pop + F1*(x2 - x3) + F2*(bm - pop);
        case 5 % DE/rand/2
            vi = x1 + F1*(x2 - x3 + x4 - x5);
        case 6 % DE/best/2
            vi = bm + F1*(x2 - x3 + x4 - x5);
    end
    
    %% Recombination
    rMat  = (rand(NP, D) < CR) | (repmat(randi(D, NP, 1),1,D) == compVec);
    nrMat = not(rMat);
    ui    = nrMat.*pop + rMat.*vi;
    ui    = ui./(2*sum(ui,2));
    ui    = round( ui.*2.^( 15 - ceil(log2(max( abs(ui), [], 2 ))) ) );
    
    %% Score these so we can take the best of each population
    [tmpscore, ~]         = TestPopulation([ui, fliplr(ui)], targInput, false); % Test the initial population
    winPop                = tmpscore <= score;      % Find the mutations that are better than original
    score(winPop)         = tmpscore(winPop);       % Replace the old score with the better scores
    pop(winPop,:)         = ui(winPop,:);           % Replace the old members with the better members
    [currBstScore, bstIx] = min(score);
    bm                    = pop(bstIx,:);
    
    %% Update the figures every so often so we know how things are going
    scrCnt                = rem(cnt-1, refScore)+1; % loop a counter
    if(scrCnt==1), bstScore = NaN(refScore, 1); allScore = NaN(refScore, NP);end
    bstScore(scrCnt)      = currBstScore;
    allScore(scrCnt,:)    = score.';
    
    if(rem(cnt, refVal)==0)
        figure(100);
        subplot(221),plot(10*log10(allScore),'linewidth',2),hold on,
        plot(10*log10(bstScore),'-m','linewidth',2),hold off,grid on%,ylim([10*log10(min(min(allScore)))-10,10+10*log10(max(max(allScore)))])
        subplot(222),plot(-diff(10*log10(allScore)),'linewidth',2),hold on,
        plot(-diff(10*log10(bstScore)),'-m','linewidth',2),hold off,grid on%,ylim([10*log10(min(min(allScore)))-10,10+10*log10(max(max(allScore)))])
        subplot(212),plot(log2(abs([pop,fliplr(pop)].')),'linewidth',2);hold on;
        plot(log2(abs([bm,fliplr(bm)])),'-m','linewidth',2);hold off;grid on
        title(sprintf('Current Population, Iter=%g', cnt));
        set(figure(100),'Position',[300 -460 560 400 ]);
        %     drawnow
        
        TestPopulation([bm, fliplr(bm)], targInput, true);
    end
    
    %% Index the loop counter
    cnt = cnt + 1;
end
