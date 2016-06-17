classdef Node < handle
    
    properties (SetAccess = private)
        CPT
        inNodes
        outNodes
        
        counterInNew
        counterOut
        
        newInSeparator
        oldInSeparator
        
        OutSeparator
    end
    
    methods        
        function node = Node(cpt, inNodes, outNodes)
            
            node.CPT = cpt;
            node.inNodes = inNodes;
            node.outNodes = outNodes;
            
            node.oldInSeparator = cell(length(inNodes),1);
            node.newInSeparator = cell(length(inNodes),1);
            node.counterInNew = 0;
            
            for i=1:length(inNodes)
                node.oldInSeparator{i} = CPT([],[]);
            end
            
            node.OutSeparator = cell(length(outNodes),1);
            node.counterOut = 0;
        end
        
        
        function cpt = getCPT(node)
            cpt = node.CPT.getCopy();
        end
        
        %function which return copy of this node, contatining only CPT
        function nodeCopy = getPartialCopy(node)
            cpt = node.getCPT();
            nodeCopy = Node(cpt, [], []);
        end
        
        
        function addNewInSeparator(node, CPT, ind)
            node.counterInNew = node.counterInNew + 1;
            node.newInSeparator{ind} = CPT;
        end
        
        
        function addOutSeparator(node, CPT, ind)
            node.counterOut = node.counterOut + 1;
            node.OutSeparator{ind} = CPT;
        end        
        
        function reverse(node)
            %switch inputs and outputs
            %and clear input and output containers from CPTs
            
            tmp = node.inNodes;
            node.inNodes = node.outNodes;
            node.outNodes = tmp;
            
            node.oldInSeparator = node.OutSeparator;            
            node.newInSeparator = cell(length(node.inNodes), 1);                                    
            
            node.OutSeparator = cell(length(node.outNodes), 1);
            
            node.counterInNew = 0;
            node.counterOut = 0;
        end        
    end
end