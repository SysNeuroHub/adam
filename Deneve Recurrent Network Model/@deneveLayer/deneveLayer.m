classdef deneveLayer < handle
    
    properties
        name@char;
        inLayer@deneveLayer;    %The handle(s) of the input layer(s)
        w@cell;                 %Cell array of weight vectors/matrices, one for each input layer, used to pool the activity of its input layer
        r@double;               %Activity of each unit       
    end
    
    properties (Dependent)
        size;
        nInLayers;
        nUnits;
    end
    
    methods
        function sz = get.size(o)
            sz = size(o.r); %#ok<CPROP>
        end
        function n = get.nInLayers(o)
            n = numel(o.inLayer);
        end
        function n = get.nUnits(o)
            n = prod(o.size);
        end
    end
    
    methods
        function o = deneveLayer(name, dim1,dim2)
            %Constructor.
            o.name = name;
            o.r = zeros(dim1,dim2);
        end
        
        function o = addInput(o,layer)
            %Add one (or more) input layers that feed this layer
            o.inLayer = horzcat(o.inLayer,layer);
        end
        
        function o = initialise(o,r)
            %Set the initial state of the layer (i.e. activity of each unit)
            o.r = r;
        end
        
        function o = update(o)
            %Pool the input from this layer's input layer(s)
            o = poolInput(o);
            
            %Perform squaring and divisive normalisation
            o = transfer(o);
        end
        
        function o = poolInput(o)     
            %Pool the input to each unit, across all input layers (o.inLayers is > 1 only for the hidden layer)
            %Uses linear indexing so it can be blind to whether this is a 1D or 2D layer
            rIn = zeros(o.size);
            for u=1:o.nUnits
                for i=1:o.nInLayers
                    thisPool = o.w{u}{i}.*o.inLayer(i).r;
                    rIn(u) = rIn(u) + sum(thisPool(:));
                end
            end
            o.r = rIn;
        end
        
        function o = transfer(o)
            S = 0.1;
            mu = 0.002;
            o.r = o.r.^2;
            o.r = o.r/(S+mu.*sum(o.r(:)));
        end
    end
end