
classdef CPT < handle
    
    properties (SetAccess = private)
        table
        dims
        varIDs
        unit
        logSpace
        logValue
    end
    
    
    methods
        function cpt = CPT(vars, dims, type)
                                    
            if isempty(vars) && isempty(dims)
                cpt.unit = 1; %unit CPT
            else                
                cpt.table = zeros(dims);
                cpt.dims = dims;
                cpt.varIDs = vars;
                cpt.unit = 0;
            end
                        
            if nargin == 2
                cpt.logSpace = 0; %regular cpt                                
            else
                cpt.logSpace = type; %cpt in log space
            end
        end
        
        
        function cpt_copy = getCopy(cpt)
            
            cpt_copy = CPT([],[]);
            
            cpt_copy.table = cpt.table;
            cpt_copy.dims = cpt.dims;
            cpt_copy.varIDs = cpt.varIDs;
            cpt_copy.unit = cpt.unit;
            cpt_copy.logSpace = cpt.logSpace;
            cpt_copy.logValue = cpt.logValue;
        end
        
        
        function fill(cpt, data)
            cpt.table(:) = data;                        
        end
        
        
        
        function multiply(cpt1, cpt2)
            
            if cpt1.unit == 1
                assert(1==0);
            elseif cpt2.unit == 1
                %do nothing
            else
                [resTable1 resTable2 resDims, resVars] = embedTables(cpt1, cpt2);
              
                %in log space
                cpt1.table = resTable1 + resTable2;
                
                cpt1.dims = resDims;
                cpt1.varIDs = resVars;
            end
        end
        
        
        function divide(cpt1, cpt2)
            
            if cpt1.unit == 1
                assert(1==0);
            elseif cpt2.unit == 1
                %do nothing
            else
                [resTable1 resTable2 resDims, resVars] = embedTables(cpt1, cpt2);
                
                %perform division of two CPTs
                
                %in log space
                ind_inf = find(resTable2 == -inf);
                resTable2(ind_inf) = 0;
                resTable1(ind_inf) = -inf;                
                cpt1.table = resTable1 - resTable2;                
                
                cpt1.dims = resDims;
                cpt1.varIDs = resVars;
            end
        end
          
        
        
        function marginalize(cpt, var)
                       
            %in log space             
            for i=1:length(var)
                varInd = find(cpt.varIDs == var(i));
                                               
                %size of the dim to be summed out
                N = size(cpt.table, varInd);
                if N > 1
                    
                    dim = cpt.dims;
                    dim(varInd) = [];                    
                    result = zeros(dim);
                    
                    for k=1:numel(result)   
                        
                        elements = zeros(N,1);
                        for t = 1:N                            
                            ind(1:ndims(cpt.table)) = {':'};
                            ind(varInd) = {t};
                            
                            tbl = cpt.table(ind{:});
                            tbl = tbl(:);
                            elements(t) = tbl(k);
                        end
                        
                        %do log of sum as a special form of sum of logs
                        result(k) = cpt.logsum(elements);
                    end                    
                    cpt.table = result;
                    
                else
                    %just remove this dimension
                    order = 1:length(cpt.varIDs);
                    order(varInd) = [];
                    order = [varInd order];
                    cpt.table = permute(cpt.table, order);
                    s = size(cpt.table);
                    cpt.table = reshape(cpt.table, [s(2:end) 1]);                    
                end
                                
                cpt.varIDs(varInd) = [];
                cpt.dims(varInd) = [];
            end
        end
        
        
        
        function maxind = maximize(cpt, var)
                    
            maxind = [];
            
            %in log space             
            for i=1:length(var)
                varInd = find(cpt.varIDs == var(i));
                                               
                %size of the dim over which we determine max
                N = size(cpt.table, varInd);
                if N > 1
                    
                    dim = cpt.dims;
                    dim(varInd) = [];                    
                    result = zeros(dim);
                    maxind = zeros(dim);
                    
                    for k=1:numel(result)   
                        
                        elements = zeros(N,1);
                        for t = 1:N                            
                            ind(1:ndims(cpt.table)) = {':'};
                            ind(varInd) = {t};
                            
                            tbl = cpt.table(ind{:});
                            tbl = tbl(:);
                            elements(t) = tbl(k);
                        end
                        
                        %do log of sum as a special form of sum of logs
                        [m, in] = max(elements);
                        result(k) = m;
                        maxind(k) = in;
                    end                    
                    cpt.table = result;
                    
                else
                    %just remove this dimension
                    order = 1:length(cpt.varIDs);
                    order(varInd) = [];
                    order = [varInd order];
                    cpt.table = permute(cpt.table, order);
                    s = size(cpt.table);
                    cpt.table = reshape(cpt.table, [s(2:end) 1]);                    
                end
                                
                cpt.varIDs(varInd) = [];
                cpt.dims(varInd) = [];
            end
        end
        
        
        
        function log_nrm = normalize(cpt)
            
            tbl = reshape(cpt.table, numel(cpt.table), 1);
            
            nrm = sum(tbl);
            
            if nrm == 0
                nrm = 1;
            end
           
            cpt.table = cpt.table/nrm;
            
            log_nrm = log(nrm);            
        end
        
        
        
        function log_nrm = lognormalize(cpt)
            
            tbl = reshape(cpt.table, numel(cpt.table), 1);
            
            %get log of sum of elements from logs of elements
            log_nrm = cpt.logsum(tbl);
            
            cpt.fill(tbl - log_nrm);
                       
        end
        
        
        function setEvidence(cpt, evidence)
            %evidence = [varID1 value1; varID2 value2; ...]
            for i = 1:size(evidence, 1)
                varInd = find(cpt.varIDs == evidence(i,1));
                evidenceInd = evidence(i,2);
                
                %create mask which has [0 0 .. 1 .. 0 0] on the dimension to be
                %squeezed and [1 1 1 .. 1 1] on all other dimensions
                vecs = cell(length(cpt.dims), 1);
                for j=1:length(cpt.dims)
                    if j == varInd
                        v = zeros(cpt.dims(j),1);
                        v(evidenceInd) = 1;
                        vecs{j} = v;
                    else
                        v = ones(cpt.dims(j),1);
                        vecs{j} = v;
                    end
                end
                
                nvecs = length(vecs);
                [grids{1:nvecs}] = ndgrid(vecs{:});
                Mask = grids{1};
                for j=2:nvecs
                    Mask = Mask .* grids{j};
                end
                
                %sum along the specific dimension to get result
                cpt.table = cpt.table.*Mask;
                cpt.table = sum(cpt.table, varInd);                
                
                cpt.dims(varInd) = 1;
            end
        end        

    end
        
    methods (Static)
                
        function result = logsum(array)
            
            if all(array==-inf)
                result = -inf;
            else
                [maxElem maxInd] = max(array);
                array(maxInd) = [];
                result = log1p(sum(exp(array - maxElem))) + maxElem;
            end
        end
    end
    
    
    
    methods (Hidden)
        function [resTable1 resTable2 resDims, resVars] = embedTables(cpt1, cpt2)
            
            %get common variables between 2 cpts
            commVars = intersect(cpt1.varIDs, cpt2.varIDs);
            
            %variables of the resulting cpt
            [resVars resVarLoc1 resVarLoc2] = union(cpt1.varIDs, cpt2.varIDs);
            
            %find sizes of these variables
            resDims = [cpt1.dims(resVarLoc1) cpt2.dims(resVarLoc2)];
            
            %since union sorts results, we need to adjust resDims
            [~, permInd] = sort([cpt1.varIDs(resVarLoc1) cpt2.varIDs(resVarLoc2)]);
            resDims = resDims(permInd);
            
            
            %permute dimensions of cpt1 to match dimensions of resTable1
            [cpt1.varIDs, order1] = sort(cpt1.varIDs);
            if length(order1)==1                
                cpt1.table = permute(cpt1.table, [1 2]);
            else
                cpt1.table = permute(cpt1.table, order1);
            end
            cpt1.dims = cpt1.dims(order1);
            
            
            [cpt2.varIDs, order2] = sort(cpt2.varIDs);
            if length(order2)==1
                cpt2.table = permute(cpt2.table, [1 2]);
            else
                cpt2.table = permute(cpt2.table, order2);
            end
            cpt2.dims = cpt2.dims(order2);
            
            %=== Embed first CPT ====
            %need to "embed" cpt1.table into a table of shape same as resulting table
            %location of dims of cpt1.table in the dims of resulting table
            [~, loc1, ~] = intersect(resVars, cpt1.varIDs);
            
            %also create complement of the above indices
            loc1_compl = 1:length(resDims);
            loc1_compl(loc1) = [];
            
            %create a template of table of shape same as resulting table
            tmpDims = resDims;
            tmpDims(loc1_compl) = 1;
            resTable1 = zeros(tmpDims);
            resTable1(:) = cpt1.table;
            
            %now tile up this matrix to match shape in "resDims"
            tmpDims = resDims;
            tmpDims(loc1) = 1;
            resTable1 = repmat(resTable1, tmpDims);
            
            
            %=== Embed second CPT ====
            %need to "embed" cpt2.table into a table of shape same as resulting table
            %location of dims of cpt2.table in the dims of resulting table
            [~, loc2, ~] = intersect(resVars, cpt2.varIDs);
            
            %also create complement of the above indices
            loc2_compl = 1:length(resDims);
            loc2_compl(loc2) = [];
            
            %create a template of table of shape same as resulting table
            tmpDims = resDims;
            tmpDims(loc2_compl) = 1;
            resTable2 = zeros(tmpDims);
            resTable2(:) = cpt2.table;
            
            %now tile up this matrix to match shape in "resDims"
            tmpDims = resDims;
            tmpDims(loc2) = 1;
            resTable2 = repmat(resTable2, tmpDims);
        end
                                
    end
end



