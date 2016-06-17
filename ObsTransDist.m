classdef ObsTransDist < handle
    
    properties
        As
    end
    
    methods
        function pdf = ObsTransDist(As)
            pdf.As = As;
        end
        
        function v = getValue(pdf, yt, yt_1)
            
            Nx = size(pdf.As,3);
            v = zeros(1, 1, Nx);
            
            %observation distribution already in log space
            for i=1:Nx
                v(i) = -1/2*( norm(yt - pdf.As(:,:,i)*yt_1)^2 );
            end
        end
        
        function update(pdf, A, index)
            pdf.As(:,:,index) = A;
        end
    end
end