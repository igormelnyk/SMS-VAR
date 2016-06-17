classdef ModeDurDist < handle
    
    properties        
        durations
    end
    
    methods
        function pdf = ModeDurDist(durations)
            pdf.durations = durations;
        end
        
        function v = getValue(pdf, dt, mt, dt_1)
            if dt_1 > 1
                if dt == dt_1 - 1
                    v = 1;
                else
                    v = 0;
                end
            else
                if isKey(pdf.durations, num2str(mt))
                    tmp = pdf.durations(num2str(mt));
                    v = poisspdf(dt, tmp.meanDur);
                else
                    v = log(0);
                end
            end
        end
        
        
        function v = getLogValue(pdf, dt, mt, dt_1)
            if dt_1 > 1
                if dt == dt_1 - 1
                    v = log(1);
                else
                    v = log(0);
                end
            else
                if isKey(pdf.durations, num2str(mt))
                    tmp = pdf.durations(num2str(mt));
                    v = -tmp.meanDur - gammaln(dt+1) + dt*log(tmp.meanDur);
                else
                    v = log(0);
                end
            end
        end
        
        
        function mu = getProb(pdf, mt, dt_1)
            
            if dt_1 == 1 %duration of prev mode expired
                
                %figure out what is the prob distrib of next state
                mu = pdf.durations(mt);
                
            else %same mode will continue with duration decremented by one
                
                mu = -1;
            end
        end                
    end
end