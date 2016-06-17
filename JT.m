
classdef JT < handle
    
    properties
        cliqueTree
        activeSet
        direction        
        root
        distrib
        logLikeModeSeq
        logLikeThisTimeSeries
        sepFor1
        sepFor2
        llIncrementHistory
        operation
        maxInd
        Nx
        maxIndT
        pathViterbi
        probForViterbi %temporal just to store prob to recover prob of Viterbi path
        probViterbiPath 
      
        histModeTransErr
        histModeDurErr
        histPhaseTransErr
        histVARErr        
    end
            
    
    methods
        
        function jt = JT(ModeData, SensorData, parameters, operation)
            
            %operation = 'max' or 'sum' in the max-product or sum-product
            %algorithms.
            %max-product is for Viterbi path
            %sum-product is for EM
            
            if nargin == 4
                [cliqueTree, activeSet, init, root] = constructPartialJT(ModeData, SensorData, parameters);
                
                jt.cliqueTree = cliqueTree;
                jt.activeSet = activeSet;
                jt.distrib = init;
                jt.root = root;
                jt.direction = 'forward';
                jt.logLikeModeSeq  = getLogLikeModeSequences(ModeData, parameters);
                jt.operation = operation;
                jt.Nx = parameters.Nx;
            else
                jt.cliqueTree = [];
                jt.activeSet = [];
                jt.direction = [];
                jt.root = [];
                jt.distrib = [];
                jt.logLikeModeSeq = [];
                jt.logLikeThisTimeSeries = [];
                jt.llIncrementHistory = [];
                jt.operation = [];
                jt.maxInd = [];
                jt.Nx = [];
                jt.maxIndT = [];
                jt.pathViterbi = [];
            end
        end
        
        
        %copy only nodes inside cliqueTree
        function jtCopy = getPartialCopy(jt)
            
            jtCopy = JT();
            jtCopy.cliqueTree = cell(length(jt.cliqueTree), 1);
            
            for i = 1:length(jt.cliqueTree)
                jtCopy.cliqueTree{i} = jt.cliqueTree{i}.getPartialCopy();
            end
        end
        
        %copy only data related to anomaly detection
        function res = getCopyWithAnomData(jt)                                                
            res.logLike = jt.llIncrementHistory;
            res.phaseErr = jt.histPhaseTransErr;    
        end
        
        %pass messages through system for EM algorithm
        function loglike = propagateMessages(jt)                        
                    
            loglike = 0;
            
            logLikeIncrements = zeros(1, length(jt.cliqueTree));
            indMax = zeros(jt.Nx, length(jt.cliqueTree)/2);
            probForViterbi = zeros(jt.Nx, jt.Nx, length(jt.cliqueTree)/2);
                       
            while(1)
                
                updatedActiveSet = zeros(length(jt.activeSet),1);
                c=1;
                
                for i = 1:length(jt.activeSet)                    
                    
                    %get current node
                    ind = jt.activeSet(i);                    
                    node = jt.cliqueTree{ind};
                    
                    %check to see if we got all incoming messages
                    %if so - proceed                    
                    if ~isempty(node.inNodes)
                        if node.counterInNew == length(node.inNodes)
                            
                            %multiply messages inside the node, k=2: multiply phase trans node by VAR node, k=1: multiply phase trans mode by prev phase trans mode 
                            for k = length(node.inNodes):-1:1
                                                               
                                %I am adding to node (xt xt_1 mt dt_1) the
                                %info from modes, which was computed
                                %earlier and stored in jt.logLikeModeSeq                                
                                if mod(node.inNodes(k),2)~=0 && strcmp(jt.direction, 'forward')
                                    assert(length(node.CPT.varIDs)==4);
                                    node.CPT.fill(node.CPT.table + sum(jt.logLikeModeSeq(:,ind/2)));
                                end

                                node.CPT.multiply(node.newInSeparator{k})
                                node.CPT.divide(node.oldInSeparator{k});
                                
                                %save product of VAR, phase trans. and mode nodes (without product of prev phase trans node) to recover probability of Viterbi path                                 
                                if mod(node.inNodes(k),2)~=0 && strcmp(jt.operation, 'max')
                                    assert(mod(ind,2) == 0);
                                    probForViterbi(:,:,ind/2) = squeeze(node.CPT.table)'; %transp is needed because after node.CPT.multiply the dimensions were sorted and we need reverse order
                                end
                                    
                            end
                            
                            %if direction is 'backward', this is final
                            %result, so we can normalize it to sum to 1
                            if strcmp(jt.direction, 'backward')
                                node.CPT.lognormalize();
                            end
                            
                        %else wait until all messages arrive
                        else
                            updatedActiveSet(c) = ind;
                            c=c+1;
                            continue;
                        end
                    end
                    
                    %get all outcoming nodes
                    outNodes = jt.cliqueTree(node.outNodes);
                    
                    for k = 1:length(outNodes)
                        
                        %find which variables are in itersection
                        sepVarIDs = intersect(outNodes{k}.CPT.varIDs, node.CPT.varIDs);                                                
                        marginalizedVarIDs = setdiff(node.CPT.varIDs, sepVarIDs);
                                                
                        separatorCPT = node.getCPT();
                        
                        %marginalize(sum out) or get max from all variables except sepVarIDs
                        if strcmp(jt.operation, 'sum')
                            separatorCPT.marginalize(marginalizedVarIDs);
                        else % operation='max' (for messages from VAR nodes to phase transitions this doesn't do anything)
                            maxind = separatorCPT.maximize(marginalizedVarIDs);
                        end
                                                
                        %compute a part of log-likelihood
                        if strcmp(jt.direction, 'forward') 
                            
                            %sepBefore = [sepBefore squeeze(separatorCPT.table)];
                            log_nrm = separatorCPT.lognormalize();
                            %sepAfter = [sepAfter squeeze(separatorCPT.table)];

                            loglike = loglike + log_nrm;
                            logLikeIncrements(ind) = log_nrm;
                            
                            %only for forward messages from phase trans to phase trans nodes                            
                            if strcmp(jt.operation, 'max') && mod(ind,2) == 0 
                                indMax(:, ind/2) = maxind';                                                                                                                 
                            end                            
                        end                                                                                                
                                                
                        %add this separator to outcoming node
                        sepInd = find(outNodes{k}.inNodes == ind);
                        outNodes{k}.addNewInSeparator(separatorCPT, sepInd);
                        
                        %add this separator to current node for later
                        %use as (future) incoming node (currently it is among one of outcoming nodes)
                        node.addOutSeparator(separatorCPT, k);
                    end
                    
                    nodesToInclude = node.outNodes;
                    for k=1:length(node.outNodes)
                        if any(jt.activeSet == node.outNodes(k))
                            nodesToInclude(nodesToInclude == node.outNodes(k)) = [];
                        end
                    end
                    
                    for ii=1:length(nodesToInclude)
                        updatedActiveSet(c) = nodesToInclude(ii);                    
                        c=c+1;
                    end
                end
                
                
                if nnz(updatedActiveSet) == 0
                    
                    if strcmp(jt.direction, 'forward') && strcmp(jt.operation, 'sum')
                        
                       %normalize result inside the last node at the end
                       log_nrm = jt.cliqueTree{ind}.CPT.lognormalize();                                                       
                       loglike = loglike + log_nrm;
                       logLikeIncrements(ind) = log_nrm;                               
                       
                    elseif strcmp(jt.direction, 'forward') && strcmp(jt.operation, 'max')
                        
                        %get the index of max element in the last node
                        %we will use this index to backtrack Viterbi path
                        %in the jt.pathViterbi matrix
                        sepVarIDs = intersect(jt.cliqueTree{ind}.CPT.varIDs, jt.cliqueTree{ind-1}.CPT.varIDs);
                        marginalizedVarIDs = setdiff(jt.cliqueTree{ind}.CPT.varIDs, sepVarIDs);
                                                
                        maxind = jt.cliqueTree{ind}.CPT.maximize(marginalizedVarIDs);                        
                        indMax(:, ind/2) = maxind';
                        
                        log_nrm = jt.cliqueTree{ind}.CPT.lognormalize();
                        loglike = loglike + log_nrm;
                        logLikeIncrements(ind) = log_nrm;
                        
                        [val, maxind] = max(squeeze(jt.cliqueTree{ind}.CPT.table));
                        jt.maxIndT = maxind;
                        val
                    end                                       
                    
                    break;
                else
                    updatedActiveSet = find(accumarray(nonzeros(updatedActiveSet)+1,1))-1;
                    jt.activeSet = updatedActiveSet;                    
                end
            end             
            
            if strcmp(jt.direction, 'forward')
                jt.llIncrementHistory = logLikeIncrements;
                jt.maxInd = indMax;
                jt.probForViterbi = probForViterbi;
            else
                %jt.sepBack = separBack;
            end
        end  
                
        
        %pass data through system for anomaly detection(after param learnt)
        function loglike = propagateData(jt, ModeData, SensorData, parameters)                        
                    
            loglike = 0;            
            logLikeIncrements = zeros(1, length(jt.cliqueTree)/2);
            
            Nx = parameters.Nx;
            Nm = parameters.Nm;
            Ny = parameters.Ny;
                                    
            histModeTransErr = zeros(1, length(jt.cliqueTree)/2);
            histModeDurErr = zeros(1, length(jt.cliqueTree)/2);
            histPhaseTransErr = zeros(1, length(jt.cliqueTree)/2);
            histVARErr = zeros(1, length(jt.cliqueTree)/2);
            
            PH = [];
            
            for ind = 1:length(jt.cliqueTree)/2
                %% ================== BEFORE OBSERVATION ==================
                
                %we use data in range 2:T, so ind=1 corresponds to data(2)
                %get previous mode duration
                if ind == 1
                    dt_1 = 1;
                else                    
                    dt_1 = ModeData(2, ind);
                end
                
                %get previous mode value
                mt_1 = ModeData(1, ind);
                
                %get previous VAR value
                yt_1 = SensorData(:, ind);                                                                  
                
                %===================== mode duration ======================
                if dt_1 == 1                    
                    %compute expectation of current duration
                    mt_ids = parameters.modeTransDist.transitions(num2str(mt_1)).transitionTo;
                    mt_probs = parameters.modeTransDist.transitions(num2str(mt_1)).probability;
                    
                    %expected mode duration in non-log scale
                    expect_dur = 0;
                    for ii=1:length(mt_ids)
                        dt = parameters.modeDurDist.durations(num2str(mt_ids(ii))).meanDur;
                        expect_dur = expect_dur + mt_probs(ii)*dt;
                    end
                else
                    %same mode is propagated with decremented counter
                    expect_dur = dt_1 - 1;
                end
                
                %histModeDurPrior(ind) = expect_dur;
                modeDurPrior = expect_dur;
                
                
                %===================== phase transition ===================                
                if dt_1 == 1
                    mt_ids = parameters.modeTransDist.transitions(num2str(mt_1)).transitionTo;
                    mt_probs = log(parameters.modeTransDist.transitions(num2str(mt_1)).probability);%=p(mt|mt-1,dt-1)
                    
                    p1 = parameters.phaseTransDist.T; %p(xt|xt-1,mt,dt-1)
                    
                    %multiply(add in log space) with p(mt|mt-1,dt-1)
                    tmp1 = zeros(Nx, Nx, length(mt_ids)); %=p(xt,mt|xt-1,dt-1)                    
                    for ii=1:length(mt_ids)
                        tmp1(:,:,ii) = p1(:,:,mt_ids(ii)) + mt_probs(ii);
                    end
                                                          
                    
                    %prepare p(xt-1|y1..yt-1, m1..mt-1, d1..dt-1)
                    if ind == 1
                        p2 = zeros(Nx,1);
                    else                        
                        node = jt.cliqueTree{2*ind};
                        p2 = squeeze(node.newInSeparator{1}.table);
                    end
                                        
                    %reshape p into form of tmp1
                    tmp2 = repmat(p2', [Nx, 1, length(mt_ids)]);
                    
                    %multiply(add in log space) p(xt,xt-1,mt|y1..yt-1, m1..mt-1, d1..dt-1)
                    tmp3 = tmp1 + tmp2;
                    
                    %sum over mt (log-sum) p(xt,xt-1|y1..yt-1, m1..mt-1, d1..dt-1)
                    sumed = zeros(Nx, Nx);
                    for ii=1:Nx
                        for jj=1:Nx
                            sumed(ii, jj) = CPT.logsum(squeeze(tmp3(ii,jj,:)));
                        end
                    end
                    
                    %sum over xt-1 (log-sum) p(xt|y1..yt-1, m1..mt-1, d1..dt-1)
                    phaseDist = zeros(Nx, 1);
                    for ii=1:Nx
                        phaseDist(ii) = CPT.logsum(sumed(ii,:));
                    end
                    
                    %normalize
                    nrm = CPT.logsum(phaseDist);            
                    pbefore = phaseDist - nrm;
                else
                    %same phase continues
                    node = jt.cliqueTree{2*ind};
                    pbefore = squeeze(node.newInSeparator{1}.table); %since p(xt|xt-1,mt,dt-1) and p(mt|mt-1,dt-1) are Delta functions
                end
                
                %pbefore = p(xt|y1..yt-1, m1..mt-1, d1..dt-1)                
                %histPhaseTransPrior(:, ind) = pbefore;
                phaseTransPrior = pbefore;
                                
                %===================== VAR observation ====================                
                As = parameters.obsTransDist.As;
                tmp = zeros(Ny,1);
                
                for ii=1:Nx
                    tmp = tmp + As(:,:,ii)*yt_1*exp(phaseDist(ii));                    
                end                    
                %histVARPrior(:, ind) = tmp;
                varPrior = tmp;
                
                %% =================== AFTER OBSERVATION ===================
                %we use data in range 2:T, so ind=1 corresponds to data(2)
                
                %get current mode duration
                dt = ModeData(2, ind+1);
                
                %get current mode value
                mt = ModeData(1, ind+1);
                
                %get current VAR value
                yt = SensorData(:, ind+1);
                
                %===================== mode transition ====================
                histModeTransErr(ind) = jt.logLikeModeSeq(2,ind);
                
                
                %===================== mode duration ======================
                histModeDurErr(ind) = abs(modeDurPrior - dt);
                
                %===================== phase transition ===================                
                if dt_1 == 1
                    
                    p1 = jt.logLikeModeSeq(2,ind); %p(mt=i|mt-1,dt-1)                    
                    phaseDist = parameters.phaseTransDist.T(:,:,mt); %p(xt|xt-1,mt=i,dt-1=1)                    
                    
                    
                    %multiply(add in log space) p(mt=i|mt-1,dt-1) with p(xt|xt-1,mt=i,dt-1=1)
                    tmp1 = phaseDist + p1;%=p(xt,mt=i|xt-1,mt-1,dt-1) 
                    
                    %prepare p(xt-1|y1..yt-1, m1..mt-1, d1..dt-1)
                    if ind == 1
                        p = zeros(Nx,1);
                    else                        
                        node = jt.cliqueTree{2*ind};
                        p = squeeze(node.newInSeparator{1}.table);
                    end
                                        
                    %reshape p into form of tmp1
                    tmp2 = repmat(p', [Nx, 1]);
                    
                    %multiply(add in log space)
                    tmp3 = tmp1 + tmp2;%=p(xt,mt=i,xt-1|y1..yt-1, m1..mt-1, d1..dt-1) 
                                        
                    %sum over xt-1 (log-sum)
                    phaseDist = zeros(Nx, 1);
                    for ii=1:Nx
                        phaseDist(ii) = CPT.logsum(tmp3(ii,:));
                    end
                    
                    %normalize
                    nrm = CPT.logsum(phaseDist);            
                    phaseDist = phaseDist - nrm;
                else
                    %same phase continues
                    node = jt.cliqueTree{2*ind};
                    phaseDist = squeeze(node.newInSeparator{1}.table); %=p(xt|y1..yt-1,m1..mt=i,d1..dt-1)
                end
                
                
                obsDist = squeeze(parameters.obsTransDist.getValue(yt, yt_1)); %p(yt=i|yt-1,xt)
                nrm = CPT.logsum(obsDist);
                obsDist = obsDist - nrm;
                
                
                %multiply p(yt|yt-1,xt) and p(xt|y1..yt-1,m1..mt=i,d1..dt-1)=phaseDist
                pafter = obsDist + phaseDist; %=p(xt,yt=i|y1..yt-1,m1..mt=i,d1..dt-1)                
                
                %normalize to get p(xt|y1..yt-1,yt=i,m1..mt=i,d1..dt-1)
                nrm = CPT.logsum(pafter);
                pafter = pafter - nrm;    
                
                pafter(pafter == -inf) = -1e+10;
                phaseTransPrior(phaseTransPrior == -inf) = -1e+10;
                      
                %kompute KL divergence btw two distributions
                res = 0;
%                 for ii=1:Nx
%                    res = res + exp(phaseTransPrior(ii))*(phaseTransPrior(ii) - pafter(ii));
%                 end                
                for ii=1:Nx
                    res = res + exp(pafter(ii))*(pafter(ii) - phaseTransPrior(ii));
                end

                                
                histPhaseTransErr(ind) = res;
                                
                %===================== VAR observation ====================                
                histVARErr(ind) = norm(varPrior - yt);
                                
                
                %% ================ PASS DATA FORWARD =====================
                
                nodeVAR = jt.cliqueTree{2*ind-1};
                nodeTrans = jt.cliqueTree{2*ind};
                                
                sepVarIDs = intersect(nodeTrans.CPT.varIDs, nodeVAR.CPT.varIDs);                                                
                marginalizedVarIDs = setdiff(nodeVAR.CPT.varIDs, sepVarIDs);
                                                
                separatorCPT = nodeVAR.getCPT();
                        
                %marginalize(sum out) all variables except sepVarIDs
                separatorCPT.marginalize(marginalizedVarIDs);
                
                %add this separator to outcoming node
                if ind == 1
                    nodeTrans.addNewInSeparator(separatorCPT, 1);
                else
                    nodeTrans.addNewInSeparator(separatorCPT, 2);
                end
                
                %multiply messages inside the nodeTrans
                
                %!This is a HACK! I am adding to node (xt xt_1 mt dt_1) the
                %info from modes, which was computed earlier and stored in jt.logLikeModeSeq
                nodeTrans.CPT.fill(nodeTrans.CPT.table + sum(jt.logLikeModeSeq(:,ind)));
                
                for k=1:length(nodeTrans.newInSeparator)
                    nodeTrans.CPT.multiply(nodeTrans.newInSeparator{k});
                end
                
                if ind < length(jt.cliqueTree)/2                
                    nodeTransNext = jt.cliqueTree{2*ind+2};
                                        
                    %find which variables are in intersection
                    sepVarIDs = intersect(nodeTransNext.CPT.varIDs, nodeTrans.CPT.varIDs);
                    marginalizedVarIDs = setdiff(nodeTrans.CPT.varIDs, sepVarIDs);
                    
                    separatorCPT = nodeTrans.getCPT();
                    
                    %marginalize(sum out) all variables except sepVarIDs
                    separatorCPT.marginalize(marginalizedVarIDs);
                    
                    %normalize
                    log_nrm = separatorCPT.lognormalize();                           
                    loglike = loglike + log_nrm;
                    logLikeIncrements(ind) = log_nrm;

                    nodeTransNext.addNewInSeparator(separatorCPT, 1);
                    
                    PH = [PH squeeze(separatorCPT.table)];                                
                else
                    %simply normalize result inside the last node
                    log_nrm = nodeTrans.CPT.lognormalize();                                                       
                    loglike = loglike + log_nrm;
                    logLikeIncrements(ind) = log_nrm;
                end
            end
            
            jt.llIncrementHistory = logLikeIncrements;
            jt.histModeTransErr = histModeTransErr;
            jt.histModeDurErr = histModeDurErr;
            jt.histPhaseTransErr = histPhaseTransErr;
            jt.histVARErr = histVARErr;
        end
        
        
        function logLike = run(jt)
            
            %send messages forward
            logLike = jt.propagateMessages();
            
            %store the log likelihood of this time series
            jt.logLikeThisTimeSeries = logLike;
            
            %reverse node direction
            for i=1:length(jt.cliqueTree)
                if isempty(jt.cliqueTree{i})
                    continue;
                end
                jt.cliqueTree{i}.reverse();
            end
            jt.direction = 'backward';            
            jt.activeSet = length(jt.cliqueTree);
                                    
            %send messages backward
            jt.propagateMessages();                                  
        end  
        
                
        function propagateForward(jt) 
            
            loglike = 0;
            
            logLikeIncrements = zeros(1, length(jt.cliqueTree));            
            separFor1 = [];
            separFor2 = [];
            nodeMult0 = [];
            nodeMult1 = [];
            nodeMult2 = [];
                                    
            while(1)
                
                updatedActiveSet = zeros(length(jt.activeSet),1);
                c=1;
                
                for i = 1:length(jt.activeSet)                    
                    
                    %get current node
                    ind = jt.activeSet(i);                    
                    node = jt.cliqueTree{ind};
                    
                    %check to see if we got all incoming messages
                    %if so - proceed                    
                    if ~isempty(node.inNodes)
                        if node.counterInNew == length(node.inNodes)
                            
                            %multiply messages inside the node, k=2: multiply phase trans node by VAR node, k=1: multiply phase trans mode by prev phase trans mode 
                            for k = 1:length(node.inNodes)
                                    
                                if k==1 
                                  nodeMult0 = cat(3, nodeMult0, squeeze(node.CPT.table));
                                end
                                
                                %---------
                                node.CPT.multiply(node.newInSeparator{k})                                
                                node.CPT.divide(node.oldInSeparator{k}); 
                                %---------
                                                                                                  
                                if k==1 && ind~=2
                                  separFor1 = [separFor1 squeeze(node.newInSeparator{k}.table)];
                                  nodeMult1 = cat(3, nodeMult1, squeeze(node.CPT.table)');
                                end                                
                                
                                if (k==1 && ind==2) || k==2 
                                  separFor2 = [separFor2 squeeze(node.newInSeparator{k}.table)];
                                  nodeMult2 = cat(3, nodeMult2, squeeze(node.CPT.table)');
                                end
                            end                                                        
                            
                        %else wait until all messages arrive
                        else
                            updatedActiveSet(c) = ind;
                            c=c+1;
                            continue;
                        end
                    end
                    
                    %get all outcoming nodes
                    outNodes = jt.cliqueTree(node.outNodes);
                    
                    for k = 1:length(outNodes)
                        
                        %find which variables are in itersection
                        sepVarIDs = intersect(outNodes{k}.CPT.varIDs, node.CPT.varIDs);                                                
                        marginalizedVarIDs = setdiff(node.CPT.varIDs, sepVarIDs);
                                                
                        separatorCPT = node.getCPT();
                        
                        %marginalize(sum out) all variables except sepVarIDs
                        separatorCPT.marginalize(marginalizedVarIDs);
                                                
                        %compute a part of log-likelihood
                        %sepBefore = [sepBefore squeeze(separatorCPT.table)];
                        log_nrm = separatorCPT.lognormalize();
                        %sepAfter = [sepAfter squeeze(separatorCPT.table)];
                        
                        loglike = loglike + log_nrm;
                        logLikeIncrements(ind) = log_nrm;
                                                
                        %add this separator to outcoming node
                        sepInd = find(outNodes{k}.inNodes == ind);
                        outNodes{k}.addNewInSeparator(separatorCPT, sepInd);
                        
                        %add this separator to current node for later
                        %use as (future) incoming node (currently it is among one of outcoming nodes)
                        node.addOutSeparator(separatorCPT, k);
                    end
                    
                    nodesToInclude = node.outNodes;
                    for k=1:length(node.outNodes)
                        if any(jt.activeSet == node.outNodes(k))
                            nodesToInclude(nodesToInclude == node.outNodes(k)) = [];
                        end
                    end
                    
                    for ii=1:length(nodesToInclude)
                        updatedActiveSet(c) = nodesToInclude(ii);                    
                        c=c+1;
                    end
                end                
                
                if nnz(updatedActiveSet) == 0                    
                       %normalize result inside the last node at the end
                       log_nrm = jt.cliqueTree{ind}.CPT.lognormalize();                                                       
                       loglike = loglike + log_nrm;
                       logLikeIncrements(ind) = log_nrm;                       
                       break;
                else
                    updatedActiveSet = find(accumarray(nonzeros(updatedActiveSet)+1,1))-1;
                    jt.activeSet = updatedActiveSet;                    
                end
            end             
            
            jt.sepFor1 = separFor1;
            jt.sepFor2 = separFor2;
            jt.llIncrementHistory = logLikeIncrements;
        end
                                        
                
        function T = getAllTrans(jt)
            
            T = cell(1,1);            
            c = 1;
            for i=1:length(jt.cliqueTree)
                if isempty(jt.cliqueTree{i}) || mod(i,2) ~= 0
                    continue;
                end
                
                tbl = squeeze(jt.cliqueTree{i}.CPT.table);
                T{c} = tbl;
                c = c+1;
            end                        
        end
        
        
        function T = getAllObs(jt)
            
            T = [];            
            for i=1:length(jt.cliqueTree)
                if isempty(jt.cliqueTree{i}) || mod(i,2) == 0
                    continue;
                end
                
                tbl = squeeze(jt.cliqueTree{i}.CPT.table);
                T = [T tbl];
            end
        end                        
                
    end            
end



%function to get loglikelihood of mode sequences
function logLike = getLogLikeModeSequences(ModeData, parameters)

T = size(ModeData, 2);

logLike = zeros(2, T-1);

for i=2:T
    
    mt_1 = ModeData(1, i-1);
    mt = ModeData(1, i);
        
    dt_1 = ModeData(2, i-1);
    dt = ModeData(2, i);
    
    logLike(1,i-1) = parameters.modeDurDist.getLogValue(dt, mt, dt_1);
    logLike(2,i-1) = log(parameters.modeTransDist.getValue(mt, mt_1, dt_1));
end
end



%a junction tree with implicit modes, for which cliques are not created
%explicitly
function [cliqueTree, activeSet, init, root] = constructPartialJT(ModeData, SensorData, parameters)

%ModeData: 2 x T x Nm
%SensorData: Ns x T

Nx = parameters.Nx;

%id of root node 
root = 2;

%get initial values of parameters
init.obsTransDist = parameters.obsTransDist;
init.modeDurDist = parameters.modeDurDist;
init.modeTransDist = parameters.modeTransDist;
init.phaseTransDist = parameters.phaseTransDist;

%length of time series
T = size(ModeData, 2);

%number of variables in a time slice in graphical model
M = 4;

cliqueTree = cell(2*(T-1), 1);

for i=2:T
    
    %=========== yt yt_1 xt ===============================================
    vars = [(i-1)*M+1, (i-2)*M+1, (i-1)*M+2];
    dims = [1, 1, Nx];
    
    %values of yt and yt_1
    yt = SensorData(:,i);
    yt_1 = SensorData(:,i-1);
    
    %cpt in log scale
    cpt = CPT(vars, dims);   
    cpt.fill(init.obsTransDist.getValue(yt, yt_1));
    
    inNodes = [];
    outNodes = (i-2)*2+2;
    
    node = Node(cpt, inNodes, outNodes);
    cliqueTree{(i-2)*2+1} = node;
    
    
    
    %=========== xt xt-1 mt dt-1 ==========================================
    vars = [(i-1)*M+2, (i-2)*M+2, (i-1)*M+3, (i-2)*M+4 ];
    dims = [Nx, Nx, 1, 1];
    
    %values of mode at time step i
    mt = ModeData(1, i); %add one for Matlab indexing
    
    %values of Nm mode durations at time step i-1
    if i==2
        %for first time step make duration=1 to force transition
        dur = 1;
    else        
        dur = ModeData(2, i-1);
    end
    
    %cpt in log scale
    cpt = CPT(vars, dims);
    
    %phase transition is already in LOG space
    if init.phaseTransDist.logSpace
        cpt.fill(init.phaseTransDist.getValue(mt, dur));
    else
        cpt.fill(log(init.phaseTransDist.getValue(mt, dur)));
    end
    
    %check that the sizes agree
    assert(all(size(zeros(cpt.dims)) == size(cpt.table)));
    
    if i==2
        inNodes = (i-2)*2+1;
        outNodes = (i-1)*2+2;
    elseif i==T
        %incoming nodes
        inNodes = [(i-3)*2+2, (i-2)*2+1];
        outNodes = [];
    else
        inNodes = [(i-3)*2+2, (i-2)*2+1]; 
        outNodes = (i-1)*2+2;
    end
    
    node = Node(cpt, inNodes, outNodes);
    cliqueTree{(i-2)*2+2} = node;
    
    %============ mt mt_1 dt_1 & dt mt dt_1 ===============================
    %these nodes will not be put in junction tree. they will be precomputed
    %once and used to multiply nodes 2, 4, 6, ... , 2*(T-1) later during
    %propagation of messages    
end

%indeces of nodes from which we will start progating messages
%set initial active nodes
activeSet = 1:2:2*(T-1);
end


%function to construct junction tree for a model, which provides all details
%(duration and transition of modes is written out explicitely for each step)
function jt = constructFullJT(ModeData, SensorData)

%ModeData: 2 x T x Nm
%SensorData: Ns x T

%id number of root node
root = 2;

%number of hidden phases
Nx = 10;

%number of modes
Nm = size(ModeData, 3);

%number of sensors
Ny = size(SensorData, 1);

%get initial values of parameters
init = initializeDist(ModeData, Nx, Nm, Ny);

%length of time series
T = size(ModeData, 2);

%indeces of nodes from which we will start progating messages
jt.activeSet = zeros(T-1+Nm,1);


%number of variables in a time slice
M = 2*Nm + 2;

jt.cliqueTree = cell(M*(T-1)-1, 1);


for i = 2:T
    
    %========= yt yt_1 xt ====================
    vars = [(i-1)*M+1, (i-2)*M+1, (i-1)*M+2];
    dims = [1, 1, Nx];
    
    %values of yt and yt_1
    yt = SensorData(:,i);
    yt_1 = SensorData(:,i-1);
    
    cpt = CPT(vars, dims);
    cpt.fill(init.obsTransDist.getValue(yt, yt_1));
    
    inNodes = [];
    outNodes = (i-2)*M+2;
    
    node = Node(cpt, inNodes, outNodes);
    jt.cliqueTree{(i-2)*M+1} = node;
    
    
    
    %======== xt xt_1 m1t m2t ... mNmt =======
    vars = [(i-1)*M+2, (i-2)*M+2, (i-1)*M+3 : 2 : i*M];
    dims = [Nx, Nx, ones(1, Nm)];
    
    %values of Nm modes at time step i
    mst = squeeze(ModeData(1, i, :))' + 1; %add one for Matlab indexing
    
    cpt = CPT(vars, dims);
    cpt.fill(init.phaseTransDist.getValue(mst));
    
    %check that the sizes agree
    assert(all(size(zeros(cpt.dims)) == size(cpt.table)));
    
    if i==2
        %incoming nodes
        inNodes = [1 3:2:M M+2];
        outNodes = [];
    elseif i==T
        inNodes = (i-2)*M+1;
        outNodes = (i-3)*M+2;
    else
        inNodes = [(i-2)*M+1, (i-1)*M+2];
        outNodes = (i-3)*M+2;
    end
    
    node = Node(cpt, inNodes, outNodes);
    jt.cliqueTree{(i-2)*M+2} = node;
    
    
    
    %======== modes ========================
    for k = 1:Nm
        
        %======== mt mt_1 dt_1 =============
        
        vars = [(i-1)*M+2+(k-1)*2+1, (i-2)*M+2+(k-1)*2+1, (i-2)*M+2+(k-1)*2+2];
        dims = [1, 1, 1];
        
        %values of mt and mt_1 and dt_1
        mt = ModeData(1, i, k) + 1;
        mt_1 = ModeData(1, i-1, k) + 1;
        
        if i==2 %force the first duration to be one to evoke transition
            dt_1 = 1;
            if mt == 1
                mt_1 = 2;
            else
                mt_1 = 1;
            end
        else
            dt_1 = ModeData(2, i-1, k);
        end
        
        cpt = CPT(vars, dims);
        cpt.fill(init.modeTransDist{k}.getValue(mt, mt_1, dt_1));
        
        if i==2
            inNodes = (i-2)*M+2+(k-1)*2+2;
            outNodes = 2;
        elseif i==T
            inNodes = [];
            outNodes = (i-3)*M+2+(k-1)*2+2;
        else
            inNodes = (i-2)*M+2+(k-1)*2+2;
            outNodes = (i-3)*M+2+(k-1)*2+2;
        end
        
        node = Node(cpt, inNodes, outNodes);
        jt.cliqueTree{(i-2)*M+2+(k-1)*2+1} = node;
        
        
        
        %======== dt mt dt_1 ==============
        if i<T
            
            vars = [(i-1)*M+2+(k-1)*2+2, (i-1)*M+2+(k-1)*2+1, (i-2)*M+2+(k-1)*2+2];
            dims = [1, 1, 1];
            
            %values of mt and mt_1 and dt_1
            dt = ModeData(2, i, k);
            mt = ModeData(1, i, k) + 1;
            
            if i==2
                dt_1 = 1;
            else
                dt_1 = ModeData(2, i-1, k);
            end
            
            cpt = CPT(vars, dims);
            cpt.fill(init.modeDurDist{k}.getValue(dt, mt, dt_1));
            
            inNodes = (i-1)*M+2+(k-1)*2+1;
            outNodes = (i-2)*M+2+(k-1)*2+1;
            
            node = Node(cpt, inNodes, outNodes);
            jt.cliqueTree{(i-2)*M+2+(k-1)*2+2} = node;
        end
    end
    
    %set initial active nodes
    jt.activeSet(i-1) = (i-2)*M+1;
    if i==T
        tmp = zeros(Nm,1);
        for k=1:Nm
            tmp(k) = (i-2)*M+2+(k-1)*2+1;
        end
        jt.activeSet(i:end) = tmp;
    end
end

jt.activeSet = sort(jt.activeSet, 'descend');
jt.direction = 'forward';
jt.root = root;
end


























