classdef deneveNet < dynamicprops
    
    properties

    end
    
    properties (SetAccess = protected)
        name@char;
        layers@deneveLayer;
    end
    
    properties (Dependent)
        nLayers;                   %Number of layers
    end
    
    methods
        %get number of input layers
        function n = get.nLayers(o)
            n = numel(o.layers);
        end
    end
    
    methods
        function o = deneveNet(name)
            %Constructor.
            o.name = name;
        end
        
        function addLayer(o,layer)
            %Add a deneveLayer to the network
            addprop(o,layer.name);
            o.(layer.name) = layer;
            
            %Store the handle
            o.layers = horzcat(o.layers,layer);
        end
        
        function o = reset(o)
            %Reset the network responses to zero between simulations
            for i=1:o.nLayers
                o.layers(i).reset();
            end
        end
        
        function o = update(o)
            %Update all of the layers
            for i=1:o.nLayers
                o.layers(i).update();
            end
        end
        
        function allocLog(o,nIter,nSims)
            %Allocate memory to log network state
            for i=1:o.nLayers
                allocLog(o.layers(i),nIter,nSims);
            end
        end
    end
end