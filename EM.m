%function to perform EM

classdef EM < handle
    
    properties
        data
        parameters
        jtrees
        logSpace %perform operations on phase transitions in log space
    end
    
    methods
        
        function em = EM(Data, numPhases, flag, trueDist)
            
            %flag is 1 means perform operations on phase transitions in log space
            em.logSpace = flag;            
            
            em.data = Data;
            
            %number of hidden phases
            Nx = numPhases;
            
            %number of modes
            Nm = 2^Data{1}.numModes;
            
            %number of sensors
            Ny = size(Data{1}.SensorData, 1);
            
            if isempty(trueDist)
                em.parameters = initializeDistrib(Data, Nx, Nm, Ny, flag);
            else
                em.parameters = trueDist;
            end

            em.parameters.Nx = Nx;
            em.parameters.Ny = Ny;
            em.parameters.Nm = Nm;
        end
        
        
        function run_parallel(em, Nproc, NemIter)
            
            %clean possible leftovers
            system('rm -rf sim/res_*.mat');
            system('rm -rf sim/data_*.mat');
            system('rm -rf sim/output*.txt');
                       
            %do EM iterations
            for i=1:NemIter                 
                fprintf('EM ITERATION %i\n', i);                
                
                logLike = em.Estep_parallel(Nproc)
                                
                em.Mstep(i); 
                
                %remove res_*.mat files before next iteration
                while(1)
                    if ~isempty(dir('sim/res_*.mat'))
                        system('rm -rf sim/res_*.mat');
                        pause(2);
                    else
                        break;
                    end
                end                
            end
            
            %clean possible leftovers
            system('rm -rf sim/res_*.mat');
            system('rm -rf sim/data_*.mat');
            system('rm -rf sim/output*.txt');
            
            %Viterbi
            %em.computeViterbiPath_parallel(Nproc);
        end
        
        
        function run(em, NemIter) 
            for i=1:NemIter    
                fprintf('EM ITERATION %i\n', i);
                
                logLike = em.Estep()
                em.Mstep(i);                                                
            end            
            em.computeViterbiPath();                                    
        end
                
        
        function logLike = Estep(em)
            
            N = length(em.data);
            
            %construct N junction trees and do forward-backward pass
            em.jtrees = cell(N,1);
            
            logLike = 0;            
            for i = 1:N
                em.jtrees{i} = JT(em.data{i}.ModeData, em.data{i}.SensorData, em.parameters, 'sum');
                logLike = logLike + em.jtrees{i}.run();
            end
        end
        
        
        function logLike = Estep_parallel(em, Nproc)
                        
            N = length(em.data);
            
            em.jtrees = cell(N,1);                        
            logLike = 0;          
            ind = reshape(1:N, N/Nproc, Nproc);
            
            par = em.parameters;                    
            
            for p = 1:Nproc 
                md = cell(length(ind(:,p)),1);
                sd = cell(length(ind(:,p)),1);
                
                for i = 1:length(ind(:,p))                    
                    md{i} = em.data{ind(i,p)}.ModeData;
                    sd{i} = em.data{ind(i,p)}.SensorData;                   
                end  
                
                indFlights = ind(:,p);
                
                %-v6 for faster save and load
                save(['sim/data_', num2str(p), '.mat'], '-v6', 'md', 'sd', 'par', 'indFlights');
            end
            
            clear md;
            clear sd;
            
            for p = 1:Nproc                                                
                  system(['matlab -nosplash -nodisplay -nodesktop -r "p=', ...
                        num2str(p), '; N=', num2str(N), '; Nproc=', num2str(Nproc),...
                        '; scriptEM" >&sim/output', num2str(p), '.txt &']);                                    
            end
            
            
            %check to see if all processes are done
            prev = 0;
            fprintf('Progress: ..0');
            while(1)                
                if length(dir('sim/res_*.mat')) == N
                    fprintf('..100 Done!\n');
                    break;
                else
                    num = length(dir('sim/res_*.mat'));
                    if (100*num)/N >=  prev + 10;
                        fprintf('..%3.1f',(100*num)/N);
                    end   
                    prev = (100*num)/N;
                    pause(10);
                end
            end
            
            %check to see if each result is ready to use
            for d = 1:N                
                
                fileName = ['sim/res_', num2str(d), '.mat'];
                
                while(1)
                    try
                        if exist(fileName)                            
                            S = load(fileName);                            
                            if ~isempty(fieldnames(S)) && ~isempty(S.jt.cliqueTree)
                               break;
                            end
                        end
                    catch exception
                        fprintf('Load file: %s, message: %s\n', fileName, exception.message);                        
                    end
                    disp(['Waiting for file ', fileName, ' to be ready for loading']);
                    pause(2);
                end                                
                
                logLike = logLike + S.logL;
            end
            
            for p = 1:Nproc                
                fileName = {['sim/data_', num2str(p), '.mat'], ['sim/output', num2str(p), '.txt']};                
                for i=1:length(fileName)
                    system(['rm -rf ', fileName{i}]);
                    while(exist(fileName{i}))
                        pause(2);
                        system(['rm -rf ', fileName{i}]);
                    end
                end                
            end                          
        end
                        
        
        function Mstep(em, i)
            
            if em.logSpace
                em.updatePhaseLogTransDist(i);
            else
                em.updatePhaseTransDist();
            end
            
            em.updateObsDist(i);           
        end
        
        
        function logLike = computeNormalizedLogLike(em, data)
            
            jtree = JT(data.ModeData, data.SensorData, em.parameters, 'sum');            
            logLike = jtree.run();
            
            %normalize logLikelihood using length of data sequence
            logLike = logLike/size(data.SensorData, 2);                        
        end
                                           
        
        function updatePhaseTransDist(em)
            
            N = length(em.data);
            
            updatedModes = [];
            
            for i=1:N
                
                modes = squeeze(em.data{i}.ModeData(1,:,:))';
                durs = squeeze(em.data{i}.ModeData(2,:,:))';
                
                %start from second time stamp because of 1st order VAR
                %t=1 is considered initialization
                modes = modes(:,2:end);
                
                Nm = em.parameters.Nm;
                base = zeros(Nm,1);
                for j=0:Nm-1
                    base(Nm-j) = 2^j;
                end
                                
                %identify changes in seq
                %seq(j) = bin2dec(sprintf('%d',modes(:,j)));
                seq = base'*modes;
                changeInd = find(abs(diff(seq))>0)+1;                
                
                %add change at ind=1 (initial transition)
                changeInd = [1 changeInd];
                changeInd = unique(changeInd);

                
                %update phase transition distribution
                for k=1:length(changeInd)
                    
                    ind = changeInd(k);                    
                    mst =  modes(:, ind)'+1; %+1 for matlab indexing
                    dur = durs(:, ind)';
                                            
                    if any(updatedModes == seq(ind))
                        
                        %we encountered this mode combination before
                        T = em.parameters.phaseTransDist.getValue(mst, dur);
                        
                        %use transpose so that rows will correspond to
                        %x_t and columns to x_t-1
                        T = T + exp(squeeze(em.jtrees{i}.cliqueTree{2*ind}.CPT.table))';                          
                        
                        em.parameters.phaseTransDist.update(mst, T);
                    else
                        
                        %we encountered this mode combination for the first time
                        %(transpose is to have x_t-1 for rows and x_t for columns)
                        T = exp(squeeze(em.jtrees{i}.cliqueTree{2*ind}.CPT.table))';
                        em.parameters.phaseTransDist.update(mst, T);
                        
                        updatedModes = [updatedModes; seq(ind)];
                    end
                end                                                
            end
            
            %after all updates are done, renormalize pdf's to sum to 1
            em.parameters.phaseTransDist.renormalize();
            
        end        
                                
        
        function updatePhaseLogTransDist(em, iter)
            
            N = length(em.data);
            Nx = em.parameters.Nx;
                        
            mapObj = containers.Map;
            mapObjJtrees = containers.Map;
            
            for i=1:N
                
                modes = em.data{i}.ModeData(1,:);
                
                %start from second time stamp because of 1st order VAR
                %t=1 is considered initialization
                modes = modes(2:end);
                                                                
                %identify changes in sequence of modes
                changeInd = find(abs(diff(modes))>0)+1;
                
                %add change at ind=1 (initial transition)
                changeInd = [1 changeInd];
                changeInd = unique(changeInd);
                
                %contruct map
                %modesID: [seqID indInSeq] 
                for k=1:length(changeInd)
                                        
                    ind = changeInd(k);
                    tmp = [];
                    if ~isKey(mapObj, num2str(modes(ind)))
                        tmp.seqID = i;
                        tmp.indInSeq = ind;
                    else
                        tmp = mapObj(num2str(modes(ind)));
                        tmp(end+1).seqID = i;
                        tmp(end).indInSeq = ind;                        
                    end
                    mapObj(num2str(modes(ind))) = tmp;                    
                    
                    
                    %contruct another map
                    %seqID: [modesID]
                    tmp = [];
                    if ~isKey(mapObjJtrees, num2str(i))
                        tmp.modesID = modes(ind);
                        tmp.indInSeq = ind;
                    else
                        tmp = mapObjJtrees(num2str(i));
                        tmp(end+1).modesID = modes(ind);
                        tmp(end).indInSeq = ind;
                    end
                    mapObjJtrees(num2str(i)) = tmp;
                    
                end                                               
            end
            
            keys = mapObj.keys;            
            Nupdates = length(keys);
            
            %prepare cell of matrices to store updates
            allUpdates = cell(Nupdates, 1);
            keysArray = zeros(Nupdates, 1);
            counters = ones(Nupdates, 1);
            
            for i=1:Nupdates
                allUpdates{i} = zeros(Nx, Nx, length(mapObj(keys{i})));
                keysArray(i) = str2double(keys{i});
            end
                        
            keysJtrees = mapObjJtrees.keys;            
            NJtrees = length(keysJtrees);
            
            for i=1:NJtrees
                fileName = ['sim/res_', keysJtrees{i}, '.mat'];
                S = load(fileName); 
                
                v = mapObjJtrees(keysJtrees{i});
                
                for j=1:length(v)
                    
                    ind = find(keysArray == v(j).modesID);
                    c = counters(ind);
                    allUpdates{ind}(:,:,c) = squeeze(S.jt.cliqueTree{2*v(j).indInSeq}.CPT.table)';                    
                    counters(ind) = counters(ind) + 1;
                end
            end
            
            
            for i=1:Nupdates
                
                T = zeros(Nx, Nx);
                
                %compute log of sum from logs 
                if size(allUpdates{i},3) == 1
                    T = allUpdates{i};
                else
                    for k=1:Nx
                        for s=1:Nx
                            T(k,s) = CPT.logsum(allUpdates{i}(k,s,:));
                        end
                    end
                end                
                
                %normalize all columns to sum to one (log of sum = 0)
                for k=1:Nx
                   col = T(:,k);
                   logSum = CPT.logsum(col);
                   T(:,k) = T(:,k) - logSum;
                end
                
                %recover mode combination of this phase transition matrix
                %any seqID and indInSeq in this iteration is good (same mode combination)
                v = mapObj(keys{i});                
                seqID = v(1).seqID;
                indInSeq = v(1).indInSeq;
                
                modes = em.data{seqID}.ModeData(1,:);
                modes = modes(2:end);
                mt =  modes(indInSeq); 
                          
                %finally, update the parameters
                em.parameters.phaseTransDist.update(mt, T);
            end                       
        end        
                
        
        function updateObsDist(em, iter)
            
            N = length(em.data);            
            Nx = em.parameters.Nx;
            Ny = em.parameters.Ny;
            
            %min duration of part of time series, used in VAR
            minPartialDur = 3;
            
            %min duration of total time series for a specific mode, used in VAR
            minTotalDur = Ny+10;
                                    
            %size of chunks to use in TSQR
            chunks = 200;
            
            for k=1:Nx 
                                
                R = [];
                
                %scale = 0;
                
                for i=1:N
                    
                    X = [];
                    Y = [];
                    
                    %recover the distribution of phase
                    Lseq = length(em.data{i}.ModeData(1,:,1))-1;
                    phase = zeros(1, Lseq);                                        
                    
                    %when there are too many flights, load jtree from disk
                    fileName = ['sim/res_', num2str(i), '.mat'];                
                    S = load(fileName); 
                    
                    for j=1:Lseq
                        %tmp = exp(squeeze(em.jtrees{i}.cliqueTree{2*j-1}.CPT.table));
                        
                        tmp = exp(squeeze(S.jt.cliqueTree{2*j-1}.CPT.table));
                        
                        phase(j) = tmp(k);
                        %phase2(:,j) = tmp;
                    end
                            
                    %detect parts of non-zero data which we will use in
                    %linear regression
                    nnz_phase = phase >= 1e-5;                    
                    ind = [1 find(abs(diff(nnz_phase))>0)+1 length(nnz_phase)+1];
                    dur = diff(ind);
                    dur(end) = dur(end) + 1;
                    
                    
                    startInd = [];
                    endInd = [];   
                    
                    %collect non-zero pieces of satisfactory durations
                    for j=1:length(ind)-1
                        if nnz_phase(ind(j)) == 1 && dur(j) >= minPartialDur
                            startInd = [startInd; ind(j)];
                            endInd = [endInd; ind(j+1)-1];                        
                        end
                    end
        
                    for j=1:length(startInd)
                        
                        %weigh data by the phase probability
                        W = ones(Ny, 1) * phase(startInd(j):endInd(j));                        
                                                
                        %extract data
                        x = em.data{i}.SensorData(:, startInd(j) : endInd(j));
                        y = em.data{i}.SensorData(:, startInd(j)+1 : endInd(j)+1);
                        
                        X = [X; (x.*W)'];
                        Y = [Y; (y.*W)'];
                    end
                    
                    
                    %estimate VAR parameters if we have enough data
                    %solve least-squares using TSQR
                    if ~isempty(Y)
                        T = size(Y,1);
                        indeces = 1:chunks:T;
                        indeces = [indeces T+1];
                        
                        for j = 1:length(indeces)-1
                            in = indeces(j):indeces(j+1)-1;
                            Xtmp = X(in, :);
                            Ytmp = Y(in, :);
                            R = triu(qr([R; [Xtmp Ytmp]]));
                            
                            %remove zeros from R
                            numRows = min([2*Ny size(R,1)]);
                            R = R(1:numRows, 1:2*Ny);
                        end
                    end
                    
                end
                
                if size(R,1) < Ny
                    disp(size(R))
                    disp('Matrix R is small');
                end
                
                
                %update observation distribution parameters
                if ~isempty(R)
                    
                    %compute inverse only if R(1:Ny, 1:Ny) is full rank
                    if rcond(R(1:Ny, 1:Ny)) > 1e-10
                                        
                        A = R(1:Ny, 1:Ny) \ R(1:Ny, Ny+1:end);
                        A = reshape(A, Ny, Ny)';
                                                
                        %update the observation distribution
                        em.parameters.obsTransDist.update(A, k);
                    else
                        disp('Matrix R is rank def. in updateObdDist()');
                    end
                end                                
            end
        end
        
        
        function evaluateSequences(em)
            
            N = length(em.data);
            
            %construct N junction trees and do forward-backward pass
            em.jtrees = cell(N,1);
            
            for i = 1:N
                em.jtrees{i} = JT(em.data{i}.ModeData, em.data{i}.SensorData, em.parameters, 'sum');
                em.jtrees{i}.propagateData(em.data{i}.ModeData, em.data{i}.SensorData, em.parameters);                
            end
        end
        
        
        function evaluateSequences_parallel(em, Nproc)
            
            N = length(em.data);
            
            em.jtrees = cell(N,1);                        
            
            ind = reshape(1:N, N/Nproc, Nproc);
            
            par = em.parameters;                    
            
            for p = 1:Nproc 
                md = cell(length(ind(:,p)),1);
                sd = cell(length(ind(:,p)),1);
                
                for i = 1:length(ind(:,p))                    
                    md{i} = em.data{ind(i,p)}.ModeData;
                    sd{i} = em.data{ind(i,p)}.SensorData;                   
                end          
                
                indFlights = ind(:,p);
                                
                %-v6 for faster save and load
                save(['sim/data_', num2str(p), '.mat'], '-v6', 'md', 'sd', 'par', 'indFlights');
             end
            
            
            for p = 1:Nproc                                
                system(['matlab -nosplash -nodisplay -nodesktop -r "p=', ...
                        num2str(p), '; N=', num2str(N), '; Nproc=', num2str(Nproc),...
                        '; scriptEvaluateSequences" >&sim/output', num2str(p), '.txt &']);                                   
            end
            
            
            %check to see if all processes are done
            while(1)
                if length(dir('sim/res_*.mat')) == N
                    break;
                else
                    pause(10);
                end
            end
            
            %assemble results back
            for d = 1:N
                
                fileName = ['sim/res_', num2str(d), '.mat'];
                
                while(1)
                    try
                        if exist(fileName)
                            
                            S = load(fileName);
                            
                            if ~isempty(fieldnames(S)) && ~isempty(S.res.logLike)
                                break;
                            end
                        end
                    catch exception
                        fprintf('Load file: %s, message: %s\n', fileName, exception.message);                        
                    end
                    disp(['Waiting for file ', fileName, ' to be ready for loading']);
                    pause(2);
                end                
            end
            
            for p = 1:Nproc
                fileName = {['sim/data_', num2str(p), '.mat'], ['sim/output', num2str(p), '.txt']};                
                for i=1:length(fileName)
                    system(['rm -rf ', fileName{i}]);
                    while(exist(fileName{i}))
                        pause(2);
                        system(['rm -rf ', fileName{i}]);
                    end
                end                
            end                          
        end
        
        
        %utility function
        function T = getT(em)
            T = em.parameters.phaseTransDist.T;
        end
        
        
        function A = getA(em)
            A = em.parameters.obsTransDist.As;
        end
        
        
        %utility function         
        function phase = getPhases(em, k)
            
            Nx = em.parameters.Nx;
            %             N = length(em.data);
            %             phases = cell(N,1);
            
            Lseq = length(em.data{k}.ModeData(1,:))-1;
            phase = zeros(Nx, Lseq);
            
            for j=1:Lseq
                phase(:, j) = squeeze(em.jtrees{k}.cliqueTree{2*j-1}.CPT.table);
            end            
        end
        
        
        %utility function
        function phTrans = getPhaseTrans(em, k)
            
            Nx = em.parameters.Nx;
            
            Lseq = length(em.data{k}.ModeData(1,:))-1;
            phTrans = zeros(Nx, Nx, Lseq);
            
            for j=1:Lseq
                phTrans(:, :, j) = (squeeze(em.jtrees{k}.cliqueTree{2*j}.CPT.table))';
            end
        end                                                     
        
        
        function res = evaluateObsDist(em, k, obsDist)
            
            %length of time series
            T = size(em.data{k}.ModeData, 2);

            res = zeros(em.parameters.Nx, T-1);
            
            if isempty(obsDist)
                obsDist = em.parameters.obsTransDist;
            end
            
            for i=2:T
                
                %values of yt and yt_1
                yt = em.data{k}.SensorData(:,i);
                yt_1 = em.data{k}.SensorData(:,i-1);
                
                v = squeeze(obsDist.getValue(yt, yt_1));
                res(:,i-1) = v;
            end
        end
        
                
        function res = getModeChanges(em, k)
            
            Nm = em.parameters.Nm;
            
            modes = squeeze(em.data{k}.ModeData(1,2:end,:))';
                        
            base = zeros(Nm,1);
            for j=0:Nm-1
                base(Nm-j) = 2^j;
            end
            
            %identify changes in seq
            seq = base'*modes;
            changeInd = find(abs(diff(seq))>0)+1;
            
            %add change at i=1 (initial transition)
            changeInd = [1 changeInd];
            changeInd = unique(changeInd);
            
            res = zeros(Nm+1, length(changeInd));
            for i=1:length(changeInd)
                res(1, i) = changeInd(i);
                res(2:end, i) = squeeze(em.data{k}.ModeData(1,changeInd(i)+1,:))+1;
            end            
        end                
    end        
end



function init = initializeDistrib(data, Nx, Nm, Ny, flag)

%flag = 1 means perform operations on phase transitions in log space

%generate same parameters
%rng(0);

N = length(data);

%==== initialize mode duration and transition distributions =====
durations = containers.Map;
transitions = containers.Map;

for k = 1:N
    seq = data{k}.ModeData(1,:);
    
    ind = [1 find(abs(diff(seq))>0)+1 ];
    dur = diff([ind length(seq)+1]);
    
    %keep track of mode durations
    for j=1:length(ind)                        
        tmp = [];
        if ~isKey(durations, num2str(seq(ind(j))))
            tmp.meanDur = dur(j);
            tmp.numEncounters = 1;
        else
            tmp = durations(num2str(seq(ind(j))));
            tmp.meanDur = tmp.meanDur + dur(j);
            tmp.numEncounters = tmp.numEncounters + 1;
        end
        durations(num2str(seq(ind(j)))) = tmp;
    end
    
    %if length(ind)==1 it means the is one deterministic self-transition
    %for next loop to accomodate this case, do
    if length(ind) == 1
        ind = [ind ind+1];
    end
    
    %keep track of mode transitions
    for j=2:length(ind)           
        tmp = [];
        if ~isKey(transitions, num2str(seq(ind(j-1))))
            tmp.transitionTo = seq(ind(j));
            tmp.probability = 1;
        else
            tmp = transitions(num2str(seq(ind(j-1))));
            in = find(tmp.transitionTo == seq(ind(j)));
            if isempty(in)
                tmp.transitionTo = [tmp.transitionTo seq(ind(j))];
                tmp.probability = [tmp.probability 1];
            else
                tmp.probability(in) = tmp.probability(in) + 1;
            end
        end
        transitions(num2str(seq(ind(j-1)))) = tmp;
    end    
end

%compute mean of observed modes durations
keys = durations.keys;
Nkeys = length(keys);
for i=1:Nkeys
    
    tmp = durations(keys{i});    
    tmp.meanDur = tmp.meanDur/tmp.numEncounters;
    durations(keys{i}) = tmp;
end
init.modeDurDist = ModeDurDist(durations);


%normalize transition matrix 
keys = transitions.keys;
Nkeys = length(keys);
for i=1:Nkeys
    
    tmp = transitions(keys{i});    
    tmp.probability = tmp.probability/sum(tmp.probability);
    transitions(keys{i}) = tmp;
end
init.modeTransDist = ModeTransDist(transitions);


%==== initialize phase transition distribution ====

%initialize as random
T = rand([Nx, Nx, Nm]);
S = sum(T,1);
T = T./repmat(S, [Nx, 1, 1]);
T(isnan(T))=0;

%USE T IN LOG SPACE for phase transition distribution
%so that we don't need to transform it to log space multiple times
if flag
    init.phaseTransDist = PhaseTransDist(log(T), 1);
else    
    init.phaseTransDist = PhaseTransDist(T, 0);
end

%==== initialize observation transition distribution ====
As = randn(Ny, Ny);
As = repmat(As, [1, 1, Nx]);
%As = randn(Ny, Ny, Nx);
init.obsTransDist = ObsTransDist(As);

end















