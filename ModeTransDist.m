classdef ModeTransDist < handle
    
    properties
        transitions
    end
    
    methods
        
        function pdf = ModeTransDist(transitions)
            pdf.transitions = transitions;
        end
        
        function v = getValue(pdf, mt, mt_1, dt_1)
            
            if dt_1 > 1
                if mt == mt_1
                    v = 1;
                else
                    v = 0;
                end
            else
                if isKey(pdf.transitions, num2str(mt_1))
                    tmp = pdf.transitions(num2str(mt_1));
                    ind = find(tmp.transitionTo == mt);
                    if ~isempty(ind)
                        v = tmp.probability(ind);
                    else
                        v = 0;
                    end
                else
                    v = 0;
                end
            end
        end
        
        
        function T = getProb(pdf, mt_1, dt_1)
            
            if dt_1 == 1 %duration of prev mode expired
                
                %figure out what is the prob distrib of next state
                T = full(pdf.transitions(:,mt_1));
                
            else %same mode will continue since duration is not expired yet
                
                T = zeros(size(pdf.transitions,1),1);
                T(mt_1) = 1;
            end
        end
        
        
    end
end