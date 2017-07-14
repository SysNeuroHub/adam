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
        
        
        %Eventually x will be based on the xlim range, not the size of the
        %matrix
        function rad = x2rad(o,x)
            %Functon to convert 1-o.nUnits scale to radians
            rad = x.*(2*pi)./o.size(2);
        end
        
        function x = rad2x(o,rad)
            %Functon to convert radians to 1-o.nUnits scale
            x = rad.*o.size(2)./(2*pi);
        end
        
        function error = err(o,est,true)
            %Function to determine absolute error
            estRad = x2rad(o,est);
            truRad = x2rad(o,true);

            errRad = atan2(sin(estRad-truRad),cos(estRad-truRad));
            error = rad2x(o,errRad);
        end
        
        function o = initialise(o,r)
            %Set the initial state of the layer (i.e. activity of each unit)
            o.resp = r;
        end
        
        function o = reset(o,dim1,dim2)
            %Reset the network responses to zero between simulations
            o.resp = zeros(dim1,dim2);
        end
        
        function o = addNoise(o)
            %add noise to the response
            for i = 1:o.nUnits
            o.resp(i) = poissrnd(o.resp(i));
            end
        end
        
        function vmweight = vmwfun(o,K,sigma,v)
            %Von mises weight function
            [i,j] = meshgrid(1:o.nUnits,1:o.nUnits);
            vmweight = K.*exp((cos((i-j).*(2*pi/o.nUnits))-1)/(sigma^2))+v;
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
        
        %Plot the response
        function plotState(o,varargin)
            %Parse input - change this so that it is one input
            %('plotArgs',{}) for line styles and widths
            %Remove time from class, put in deneve
            p = inputParser;
            p.addParameter('linestyle',[]);
            p.addParameter('linewidth',[]);
            p.parse(varargin{:});
            
            %Extract values from parser
            style = p.Results.linestyle;
            width = p.Results.linewidth;
            
            %Plot response
            if isvector(o.resp);     %If 1D
                plot(o.resp,style,'linewidth',width);
                hold on
                estPos = pointEstimate(o);
                maxVal = max(o.resp);
            
                plot([estPos estPos],[0 maxVal],style,'linewidth',width);
                ylim([0,100]);
            
            else                      %If 2D
                surf(o.resp); zlim([0,10])
            end
            
        end
                    
    end
end