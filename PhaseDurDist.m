classdef PhaseDurDist < handle
    
    properties        
        muArray
    end
    
    methods
        function pdf = PhaseDurDist(muArray)
            pdf.muArray = muArray;
        end
        
        function v = getValue(pdf, dt, xt, dt_1)
            if dt_1 > 1
                if dt == dt_1 - 1
                    v = 1;
                else
                    v = 0;
                end
            else
%                 v = (pdf.muArray(mt)^dt)/(factorial(dt)) * exp(-pdf.muArray(mt));
                v = poisspdf(dt, pdf.muArray(xt));
            end
        end                        
    end
end