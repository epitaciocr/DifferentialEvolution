%%
clear variables;
format long g;
tic;

%% Filter Optimization Configurations
filtnum = 5; % Filter we are trying to match
initfilt = false;
refScore = 201;
normPop = true; % Normalize and make fixed-point the members
npf = 10; % # of pops
refresh = 50; % # of iterations to output plots & text

n = 2;
% n = 20;
% targInput.ScrFnc = @(x) (mean(abs(x).^n))^(1/n);
targInput.ScrFnc = @(x) (mean(abs(x).^n));
targInput.rip = 0.05; % pband peak variation in dB
targInput.rej = 60; % dB rejection
targInput.fpass = 15;%5.2; % Passband corner freq in MHz
targInput.fstop = targInput.fpass*1.85; % targInput.fpass*9.175/4.865; % Stopband freq in MHz
% targInput.fstop = targInput.fpass*1.88; % targInput.fpass*9.175/4.865; % Stopband freq in MHz
% targInput.fstop = targInput.fpass*25.83/15.625; % Stopband freq in MHz

targInput.pbscale = 1;%1;%1e6;
targInput.sbscale = 1;%1e2;
targInput.isScl = 1; % in spec scaling
targInput.osScl = 1; % out of spec scaling

targInput.fs = 100; % reference sample frequency MHz
targInput.nTap = 64; % number of filter coeffs
targInput.fftlen = 2^10; % fft length for scoring
targInput.sbramp = [1,1]; % stop-band ramp function

% Quantization of fpass and fstop to the fft resolution
targInput.fpass = ceil( targInput.fpass/(targInput.fs/targInput.fftlen/2) )*(targInput.fs/targInput.fftlen/2);
targInput.fstop = floor( targInput.fstop/(targInput.fs/targInput.fftlen/2) )*(targInput.fs/targInput.fftlen/2);

%% Optimization Parameter Controls.
VTR = 1.e-6;            % VTR		"Value To Reach" (stop when ofunc < VTR)
XVmin = -1.01; %*ones(1,D);    % XVmin,XVmax   vector of lower and bounds of initial population
XVmax = 1.01; %*ones(1,D);
itermax = 20e6; % itermax       maximum number of iterations (generations)
F = 0.8;%0.8; % F             DE-stepsize F ex [0, 2]
CR = 0.8;  % CR            crossover probabililty constant ex [0, 1]
strategy = 1;
%        1 --> DE/best/1/exp           6 --> DE/best/1/bin
%        2 --> DE/rand/1/exp           7 --> DE/rand/1/bin
%        3 --> DE/rand-to-best/1/exp   8 --> DE/rand-to-best/1/bin
%        4 --> DE/best/2/exp           9 --> DE/best/2/bin
%        5 --> DE/rand/2/exp           else  DE/rand/2/bin

b = [1,1];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START FCN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[targInput.h_targ, ~] = freqz(b, 1, targInput.fftlen, targInput.fs);
if initfilt
    initmem = b(1:32); % Placeholder in case we want to add initialization of our DE guess
    XVmax = sqrt(mean(abs(initmem).^2));XVmin = -XVmax;
else
    initmem = 0; % Placeholder in case we want to add initialization of our DE guess
end
nTap = targInput.nTap;
if 2*round(nTap/2)==nTap
    D = nTap/2;         % D		number of parameters of the objective function
else
    error('Use an even length filter')
end

NP = 30;%npf*D; % NP            number of population members
pop = zeros(NP,D); %initialize pop to gain speed

%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

for i=1:NP
    
    pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin) + initmem;
    
    if(normPop)
        pop(i,:) = pop(i,:)./(2*sum(pop(i,:)));
        pop(i,:) = round( pop(i,:).*2.^( 15 - ceil(log2(max( abs( pop(i,:) )))) ) );
    end
end
if(initfilt)
    pop(1,:) = initmem;
    pop(1,:) = round( pop(1,:).*2.^( 15 - ceil(log2(max( abs( pop(1,:) )))) ) );
end

popold    = zeros(size(pop));     % toggle population
val       = zeros(1,NP);          % create and reset the "cost array"
bestmem   = zeros(1,D);           % best population member ever
bestmemit = zeros(1,D);           % best population member in iteration
nfeval    = 0;                    % number of function evaluations

%% -------------Evaluate the best member after initialization--------------
ibest   = 1;                      % start with first population member
[ val(1) ] = test_fir_pop_member( pop(ibest,:), targInput );
bestval = val(1);                 % best objective function value so far
nfeval  = nfeval + 1;

for i=2:NP                        % check the remaining members
    [ val(i) ] = test_fir_pop_member( pop(i,:), targInput );
    nfeval  = nfeval + 1;
    if (val(i) < bestval)           % if member is better
        ibest   = i;                 % save its location
        bestval = val(i);
    end
end
bestmemit = pop(ibest,:);         % best member of current iteration
bestvalit = bestval;              % best value of current iteration

bestmem = bestmemit;              % best member ever

%% ------DE-Minimization--------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

pm1 = zeros(NP,D);              % initialize population matrix 1
pm2 = zeros(NP,D);              % initialize population matrix 2
pm3 = zeros(NP,D);              % initialize population matrix 3
pm4 = zeros(NP,D);              % initialize population matrix 4
pm5 = zeros(NP,D);              % initialize population matrix 5
bm  = zeros(NP,D);              % initialize bestmember  matrix
ui  = zeros(NP,D);              % intermediate population of perturbed vectors
mui = zeros(NP,D);              % mask for intermediate population
mpo = zeros(NP,D);              % mask for old population
rot = (0:1:NP-1);               % rotating index array (size NP)
rotd= (0:1:D-1);                % rotating index array (size D)
rt  = zeros(NP);                % another rotating index array
rtd = zeros(D);                 % rotating index array for exponential crossover
a1  = zeros(NP);                % index array
a2  = zeros(NP);                % index array
a3  = zeros(NP);                % index array
a4  = zeros(NP);                % index array
a5  = zeros(NP);                % index array
ind = zeros(4);

iter = 1;
while ((iter < itermax) && (bestval > VTR))
    popold = pop;                   % save the old population
    
    ind = randperm(4);              % index pointer array
    
    a1  = randperm(NP);             % shuffle locations of vectors
    rt = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt = rem(rot+ind(2),NP);
    a3  = a2(rt+1);
    rt = rem(rot+ind(3),NP);
    a4  = a3(rt+1);
    rt = rem(rot+ind(4),NP);
    a5  = a4(rt+1);
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    for i=1:NP                      % population filled with the best member
        bm(i,:) = bestmemit;          % of the last iteration
    end
    
    mui = rand(NP,D) < CR;          % all random numbers < CR are 1, 0 otherwise
    
    if (strategy > 5)
        st = strategy-5;		  % binomial crossover
    else
        st = strategy;		  % exponential crossover
        mui=sort(mui');	          % transpose, collect 1's in each column
        for i=1:NP
            n=floor(rand*D);
            if n > 0
                rtd = rem(rotd+n,D);
                mui(:,i) = mui(rtd+1,i); %rotate column i by n
            end
        end
        mui = mui';			  % transpose back
    end
    mpo = mui < 0.5;                % inverse mask to mui
    
    if (st == 1)                      % DE/best/1
        ui = bm + F*(pm1 - pm2);        % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 2)                  % DE/rand/1
        ui = pm3 + F*(pm1 - pm2);       % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 3)                  % DE/rand-to-best/1
        ui = popold + F*(bm-popold) + F*(pm1 - pm2);
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 4)                  % DE/best/2
        ui = bm + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
        ui = popold.*mpo + ui.*mui;           % crossover
    elseif (st == 5)                  % DE/rand/2
        ui = pm5 + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
        ui = popold.*mpo + ui.*mui;            % crossover
    end
    
    %-----Select which vectors are allowed to enter the new population------------
    for i=1:NP
        
        % Renormalize before scoring
        if(normPop)
            ui(i,:) = ui(i,:)./(2*sum(ui(i,:)));
            ui(i,:) = round( ui(i,:).*2.^( 15 - ceil(log2(max( abs( ui(i,:) )))) ) );
        end
        
        [ tempval ] = test_fir_pop_member( ui(i,:), targInput);
        nfeval  = nfeval + 1;
        if (tempval <= val(i))  % if competitor is better than value in "cost array"
            pop(i,:) = ui(i,:);  % replace old vector with new one (for new iteration)
            val(i)   = tempval;  % save value in "cost array"
            
            %----we update bestval only in case of success to save time-----------
            if (tempval < bestval)     % if competitor better than the best one ever
                bestval = tempval;      % new best value
                bestmem = ui(i,:);      % new best parameter vector ever
            end
        end
    end %---end for imember=1:NP
    
    % Store off ever N best scores
    scrCnt = mod(iter,refScore)+1; % loop a counter
    if(scrCnt==1),bstScore=NaN(1,refScore);allScore=NaN(refScore,NP);end
    bstScore(scrCnt)=bestval;
    allScore(scrCnt,:)=val.';
    
    bestmemit = bestmem;       % freeze the best member of this iteration for the coming
    % iteration. This is needed for some of the strategies.
    
    %----Output section----------------------------------------------------------
    
    if (refresh > 0)
        if (rem(iter,refresh) == 0)
            figure(100);
            subplot(221),plot(10*log10(allScore),'linewidth',2),hold on,
            plot(10*log10(bstScore),'-m','linewidth',2),hold off,grid on%,ylim([10*log10(min(min(allScore)))-10,10+10*log10(max(max(allScore)))])
            subplot(222),plot(-diff(10*log10(allScore)),'linewidth',2),hold on,
            plot(-diff(10*log10(bstScore)),'-m','linewidth',2),hold off,grid on%,ylim([10*log10(min(min(allScore)))-10,10+10*log10(max(max(allScore)))])
            subplot(212),plot(log2(abs([pop,fliplr(pop)].')),'linewidth',2);hold on;
            plot(log2(abs([bestmem,fliplr(bestmem)])),'-m','linewidth',2);hold off;grid on
            
            subplot(212),title(sprintf('Current Population, Iter=%g', iter));
%            set(gcf,'Position',[10 550 600 400 ]);
            %
            [ ~ ] = test_fir_pop_member( bestmem, targInput, 1 );
%            set(gcf,'Position',[10 50 600 400 ]);
        end
    end
    
    iter = iter + 1;
end %---end while ((iter < itermax) ...

fprintf('Done. \n');
toc;



