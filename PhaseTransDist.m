classdef PhaseTransDist < handle
    
    properties
        T
        logSpace
    end
    
    methods
        
        function pdf = PhaseTransDist(T, logSpace)
            pdf.T = T;
            pdf.logSpace = logSpace;
        end
        
        function v = getValue(pdf, mst, durs)
                        
            Nx = size(pdf.T,1);
            
            %if any of modes changes
            if any(durs == 1)
                
                %for specific mst, select all xt, xt_1                                
                v = zeros(Nx, Nx);
                for ii=1:Nx
                    for jj=1:Nx
                        
                        ind = [ii, jj, mst];
                        
                        siz = size(pdf.T);
                        k = [1 cumprod(siz(1:end-1))];
                        ndx = 1;
                        
                        for i = 1:length(siz),
                            ndx = ndx + (ind(i)-1)*k(i);
                        end
                        
                        v(ii, jj) = pdf.T(ndx);
                    end
                end
                
            %if the same modes continue, the phase stays the same
            else
                if pdf.logSpace
                    v = log(eye(Nx));
                else                    
                    v = eye(Nx);
                end
            end                                        
        end  
        
                
        function update(pdf, mst, table)
            
            %for specific mst, set all xt, xt_1 from table
            Nx = size(pdf.T,1);
            
            for ii=1:Nx
                for jj=1:Nx
                    
                    ind = [ii, jj, mst];
                    
                    siz = size(pdf.T);
                    k = [1 cumprod(siz(1:end-1))];
                    ndx = 1;
                    
                    for i = 1:length(siz),
                        ndx = ndx + (ind(i)-1)*k(i);
                    end
                    
                    pdf.T(ndx) = table(ii, jj);
                end
            end
        end
        
        
        function renormalize(pdf)
                        
            Nx = size(pdf.T, 1);
            
            TT = pdf.T;                        
            S = sum(TT, 1);
                        
            TT = TT./repmat(S, [Nx, 1, 1]);
            %replace any NaN with non-zeros (NaN were obtained from 0/0)
            TT(isnan(TT))=0;            
            
            pdf.T = TT;            
        end
        
        
    end
end







