classdef deneveInput
    %Small class just to protect the input struct used by the deneveLayer class
    properties
        enabled@logical;
    end
    
    properties (SetAccess = protected)
        name@char;
        layer@deneveLayer;
        weights@double;
    end
    
    methods
        function o = deneveInput(layer,weights,enabled)
            if nargin < 2
                enabled = true;
            end
            
            %Constructor.
            o.name = layer.name;
            o.layer = layer;
            o.weights = weights;
            o.enabled = enabled;
        end
   
        function ind = name2ind(o,name)
           %For arrays of input objects, returns the index of the input with the specified name.
            %Used to give containers.Map-like behaviour.
           ind = find(arrayfun(@(n) strcmp(n.name,name),o));
        end
    end
end