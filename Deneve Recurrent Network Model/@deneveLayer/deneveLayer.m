classdef deneveLayer < handle
    
    properties
    %for ron    
    end
    
    properties (SetAccess = protected)
        name@char;
        input@struct;              %o.input will be a structure, with a handle, a set of weights, and a switch
        resp@double;               %Activity of each unit (response of each neuron) in the population at a given time point
        
    end
    
    properties (Dependent)
        size;                      %Size of response matrix
        nInputs;                   %Number of input layers
        nUnits;                    %Number of units in the response matrix
    end
    
    methods
        %get size of response matrix
        function sz = get.size(o)
            sz = size(o.resp);
        end
        %get number of input layers
        function n = get.nInputs(o)
            n = numel(o.input);
        end
        %get number of units in the response matrix
        function n = get.nUnits(o)
            n = prod(o.size);
        end
    end
    
    methods
        function o = deneveLayer(name, dim1,dim2)
            %Constructor.
            o.name = name;
            o.resp = zeros(dim1,dim2);
        end
        
        function o = setInput(o,layer,weights,varargin)
            %Help here.
            %This is to ensure that the input is in the correct format, and
            %allow for multiple fields after the layers and weights
            p=inputParser;
            p.addRequired('layer',@iscell);
            p.addRequired('weights',@iscell);
            p.addParameter('enabled',[], @iscell);
            p.parse(layer,weights,varargin{:});
            
            %Extract fields from input parser
            layer = p.Results.layer;
            weights = p.Results.weights;
            enabled = p.Results.enabled;
            
            %Fill the enabled field if left blank, defaults to enabled
            if isempty(enabled)
                enabled = true(1,numel(layer));
            end
            
            
            %Fill values into input structure
            nln=numel(layer);
            for i=1:nln;
                o.input(i).layer=layer{i};
                o.input(i).weights=weights{i};
                o.input(i).enabled=enabled(i);
            end
        end
        
        function o = initialise(o,r)
            %Set the initial state of the layer (i.e. activity of each unit)
            o.resp = r;
        end
        
        function o = reset(o,dim1,dim2)
            %Reset the network responses to zero between simulations
            o.resp = zeros(dim1,dim2);
        end
        
        function o = addNoise(o,n)
            %add noise to the response
            for i = 1:n
            o.resp(i) = poissrnd(o.resp(i));
            end
        end
        
        function setEnabled(o,inputName,enabled)
            %Allows enabling and disableing of input layers
            thisInput = arrayfun(@(x) strcmpi(x.layer.name,inputName),o.input);
            o.input(thisInput).enabled = enabled;
        end
        
        function o = update(o,normalise)
            %Can switch of normalisation, defaults to on
            if nargin < 2
                normalise = true;
            end
            
            %Pool the input from this layer's input layer(s)
            o = poolInput(o);
            
            %Perform squaring and divisive normalisation
            if normalise
                o = transfer(o);
            end
        end
        
        function o = poolInput(o)
            %Pool the input to each unit, across all input layers (o.inLayers is > 1 only for the hidden layer)
            %Uses linear indexing so it can be blind to whether this is a 1D or 2D layer
            rIn = zeros(o.size);
            
            for i=1:o.nInputs
                if o.input(i).enabled
                    %Only includes input if enabled
                    w = o.input(i).weights;
                    r = o.input(i).layer.resp;
                    inDims = ndims(r)-isvector(r);
                    sz = size(w);
                    %Reshape to a 2D matrix of weights by unit
                    w = reshape(w,prod(sz(1:inDims)),[]);
                    
                    for u=1:o.nUnits
                        %Sum response for a given 'neuron'
                        rIn(u) = rIn(u) + w(:,u)'*r(:);
                    end
                end
            end
            o.resp = rIn;
        end
        
        function o = transfer(o)
            %Function used to normalise responses
            S = 0.1;
            mu = 0.002;
            o.resp = o.resp.^2;
            o.resp = o.resp/(S+mu.*sum(o.resp(:)));
        end
           function estimate = pointEstimate(o)
            %Estimate position based on peak of hill
            %Preallocate matrix
            estimator = zeros(1,o.nUnits);
            %Generate matrix to sum
            for j = 1:o.nUnits
                estimator(j) = o.resp(j).*exp(1i.*((2*pi)/o.nUnits).*j);
            end
                      
            %Calculate estimate
            e = sum(estimator);
            estimate = mod(phase(e)/(2*pi)*o.nUnits,o.nUnits);
        end   
    end
end