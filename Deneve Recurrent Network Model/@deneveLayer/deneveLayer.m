classdef deneveLayer < handle
    % A population code (layer) as used in Deneve, Latham, & Pouget (2001)
    % "Efficient computation and cuintegration with noisy population codes."
    
    properties
        plotSetts = struct('lineColor','r');    %Plot settings to be used for this layer.
        xferPrms = struct('s',0.1,'mu',0.002);        %Parameters of the transfer/activation function (same for all units)
        normalise@logical = true;
        noisy@logical = true;
    end
    
    properties (SetAccess = protected)
        name@char;                              %The name of this layer
        size;                                   %Size of layer (i.e. number of units on each dimension)
        nUnits;                                 %Number of units in total
        nDims;                                  %Number of dimensions
        resp@double;                            %Activity of each unit (response of each neuron) in the population at a given time point
        input@deneveInput;                      %Struct with a layer handle, a set of pooling weights, and an enabled switch
        log@cell = {};                          %Record of previous logged layer states.
        curLog@double = 0;                      %Admin for the log.
        curPag@double = 0;                      %Ditto.
    end
    
    properties (Dependent)
        nInputs;                                %Number of input layers
    end
    
    methods
        function n = get.nInputs(o)
            n = numel(o.input);
        end
        function n = get.nUnits(o)
            n = prod(o.size);
        end 
        function n = get.nDims(o)
            n = ndims(o.resp)-isvector(o.resp);
        end
    end
    
    methods
        function o = deneveLayer(name, dim1,dim2)
            %Constructor.
            o.name = name;
            o.size = [dim1,dim2];
            o.nUnits = prod(o.size);
            o.resp = zeros(o.size);
        end
        
        function o = setInput(o,layer,varargin)
            %Add a new input to the layer. 
            %
            %INput arguments:
            %
            %   layer (required): a deneveLayer object or a [1 x N] cell array of layer objects
            %
            %Param-value pairs:
            %
            %   'weights': pooling weights for each of the input layers. A matrix or cell array.
            %   'enabled': 1 x N logical vector of length N specifying which input layers should be active
            %   'prms': 1 x N vector of parameters for the von-mises tuning function (only used if weights=='vonMises')
            
            wtValidate = @(w) (iscell(w) && (numel(w) == numel(layer)) && all(cellfun(@(x) isnumeric(x) | (ischar(x) & strcmpi(x,'vonMises')),w))) || isnumeric(w);
            p=inputParser;
            p.addRequired('layer',@(x) iscell(x) | isa(x,'deneveLayer'));
            p.addParameter('weights',repmat({'vonMises'},1,numel(layer)),wtValidate);
            p.addParameter('enabled',true(1,numel(layer)), @logical);
            p.addParameter('prms',repmat(deneveLayer.defaultVMprms('NET2NET'),1,numel(layer)));
            p.parse(layer,varargin{:});
            p = p.Results;

            %Allow non-cell input format, but change here
            if ~iscell(layer)
                layer= {layer};
            end
            if ~iscell(p.weights)
                p.weights = {p.weights};
            end

            %Add each input layer to the current layer
            nIn=numel(layer);
            %wtFuns = {@vm1Dto1D, @vm2Dto1D;@vm1Dto2D,[]}; %dim 1: my dimensionality, dim 2: input dimensionality
            wtFuns = {@vm1Dweights,@vm2Dweights};
            for i=1:nIn
                %If no weights are specified, use default weights (von-mises)
                if strcmpi(p.weights{i},'vonMises')
                    
                   %Which weight function should we use?
                   fun = wtFuns{layer{i}.nDims};
                   
                   %Use default parameters for the von-mises, but use the requested dimension
                   p.weights{i} = fun(o,layer{i},p.prms(i));
                end
                
                %Add the input object to this layer
                o.input(i) = deneveInput(layer{i},p.weights{i},p.enabled(i));
            end
        end
        
        %Eventually x will be based on the xlim range, not the size of the matrix
        function rad = x2rad(o,x)
            %Functon to convert 1-o.nUnits scale to radians
            rad = x.*(2*pi)./o.size(2);
        end
        
        function x = rad2x(o,rad)
            %Functon to convert radians to 1-o.nUnits scale
            x = rad.*o.size(2)./(2*pi);
        end
        
        function error = err(o,est,trueLoc)
            %Function to determine absolute error
            estRad = x2rad(o,est);
            truRad = x2rad(o,trueLoc);
            errRad = atan2(sin(estRad-truRad),cos(estRad-truRad));
            error = rad2x(o,errRad);
        end
        
        function setResp(o,r)
            %Set the activity of each unit to the values specified in r 
            if ~isequal(size(r),o.size) %#ok<CPROPLC>
                error('The size of the input response does not match the size of the layer.');
            end
            o.resp = r;
        end
        
        function o = reset(o)
            %Reset the network responses to zero between simulations
            o.resp = zeros(o.size);
        end
        
        function o = addNoise(o)
            %add noise to the response
            for i = 1:o.nUnits
                o.resp(i) = poissrnd(o.resp(i));
            end
        end
        
        function w = vm1Dweights(o,inLayer,prms)
            if o.nDims == 1
                w = vm1Dto1D(o,inLayer,prms);
            elseif o.nDims ==2
                w = vm1Dto2D(inLayer,o,prms);
            end
        end
        
        function w = vm2Dweights(o,inLayer,prms)
            if o.nDims == 1
                w = vm2Dto1D(inLayer,o,prms);
            elseif o.nDims ==2
                error('Not yet supported');
            end
        end
    
        function w = vm1Dto1D(o1,o2,prms)
            %Von mises weight function for 1D layer to 1D layer
            if ~isequal(o1.size,o2.size) || ~isvector(o1.size)
                error('This function only works for two 1D arrays of the same dimensionaility');
            end
            [i,j] = meshgrid(1:o1.nUnits,1:o2.nUnits);
            w = prms.k.*exp((cos((i-j).*(2*pi/o1.nUnits))-1)/(prms.sigma^2))+prms.v;
        end
        
        function w = vm1Dto2D(o1D,o2D,prms)
            
            %Von mises weight function for 1D layer to 2D layer
            %Only one dimension of the 2D layer has tuning for the 1D
            %input, specified by dim of 1, 2, or -1 for tuning along the
            %diagonal.
            [j,l,m] = meshgrid(1:o1D.nUnits,1:o2D.size(1),1:o2D.size(2));

            switch prms.dim
                case 1
                    inds = {j,l};
                case 2
                    inds = {m,l};
                case -1
                    %Diagonal requested
                    inds = {l,j+m+o1D.nUnits/2};
            end
            
            w = prms.k.*exp((cos((inds{1}-inds{2}).*(2*pi/o1D.nUnits))-1)/(prms.sigma^2))+prms.v;
        end
        
        function w = vm2Dto1D(o2D,o1D,prms)
            %Von mises weight function for 1D layer to 2D layer
            %Only one dimension of the 2D layer has tuning for the 1D
            %input, specified by dim of 1, 2, or -1 for tuning along the
            %diagonal.
            [j,l,m] = meshgrid(1:o1D.nUnits,1:o2D.size(1),1:o2D.size(2));

            switch prms.dim
                case 1
                    inds = {m,l};
                case 2
                    inds = {m,j};
                case -1
                    %Diagonal requested
                    inds = {m,j+l+o1D.nUnits/2};
            end
            
            w = prms.k.*exp((cos((inds{1}-inds{2}).*(2*pi/o1D.nUnits))-1)/(prms.sigma^2))+prms.v;
        end
        
        function allocLog(o,nPages,nLogs)
            o.log = cell(nLogs,1);
            [o.log{:}] = deal(squeeze(nan([nPages, o.size])));
        end
        
        function o = gainfun(o,trough,K,sigma)
            %For gain function K = 0.8, sigma = 0.4 and v = 0
            a = (1:o.nUnits);
            gain = 1 - (K.*exp((cos((a-trough).*(2*pi/o.nUnits))-1)./(sigma^2)));
            o.resp = o.resp.*gain;
        end
        
        function o = logState(o,varargin)
            p = inputParser;
            p.addParameter('newLog',false);
            p.parse(varargin{:});
            
            if isempty(o.log)
                error('You must prenfkdsljhf');
            end
            %Extract values from parser
            newLog = p.Results.newLog;
            
            %Use current log, or add one for new log
            if newLog
                o.curLog = o.curLog + 1;
                o.curPag = 0;
            end
            
            if o.curLog == 0
                o.curLog = 1;
            end
            
            o.curPag = o.curPag + 1;
            o.log{o.curLog}(o.curPag,:) = o.resp(:);
        end
        
        function setEnabled(o,inputName,enabled)
            %Allows enabling and disabling of input layers
            ind = name2ind(o.input,inputName);
            o.input(ind).enabled = enabled;
        end
        
        function o = update(o)

            %Pool the input from this layer's input layer(s)
            poolInput(o);
            
            %Perform squaring and divisive normalisation
            if o.normalise
                transfer(o);
            end
            
            %Use the current values as Poisson means, draw random spike counts for each neuron
            if o.noisy
                addNoise(o);
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
                    nDims = o.input(i).layer.nDims; %#ok<PROP>
                    sz = size(w); %#ok<CPROP>
                    %Reshape to a 2D matrix of weights by unit
                    w = reshape(w,prod(sz(1:nDims)),[]); %#ok<PROP>
                    r = r(:);
                    for u=1:o.nUnits
                        %Sum response for a given 'neuron'
                        rIn(u) = rIn(u) + w(:,u)'*r;
                    end
                end
            end
            o.resp = rIn;
        end
        
        function o = transfer(o)
            %Function used to normalise responses
            o.resp = o.resp.^2;
            o.resp = o.resp/(o.xferPrms.s+o.xferPrms.mu.*sum(o.resp(:)));
        end
        
        function estimate = pointEstimate(o)
            %Estimate position based on peak of hill
            %Generate matrix to sum
            estimator = zeros(1,o.nUnits);
            for j = 1:o.nUnits
                estimator(j) = o.resp(j).*exp(1i.*((2*pi)/o.nUnits).*j);
            end
            
            %Calculate estimate
            e = sum(estimator);
            estimate = mod(angle(e)/(2*pi)*o.nUnits,o.nUnits);
        end
        
        %Plot the response
        function plotState(o,varargin)
            %Parse input - change this so that it is one input
            %('plotArgs',{}) for line styles and widths
            %Remove time from class, put in deneve
            p = inputParser;
            p.addParameter('plotArgs',{}); %Param-value pairs passed to plot
            p.addParameter('ax',gca);
            p.parse(varargin{:});
            
            %Plot response
            if isvector(o.resp)     %If 1D
                plot(p.Results.ax,1:o.nUnits,o.resp,'o-','lineWidth',3,'color',o.plotSetts.lineColor,p.Results.plotArgs{:}); hold on
                estPos = pointEstimate(o);
                plot(p.Results.ax,[estPos estPos],[0 max(o.resp(:))],'lineWidth',2,'Color',o.plotSetts.lineColor);
                xlim([1,o.nUnits]);
            else                      %If 2D
                surf(o.resp);
            end
        end
        
        function heatmapMean(o,err,clim) %not sure if this will work
            merr = x2rad(o,err);
            radmean = circ_mean(merr,[],3);
            mean = rad2x(o,radmean);
            this = strcat(o.name, 'mean error');
            imagesc(mean,clim);
            title(this);
        end
        
        function heatmapStd(o,err,clim)
            merr = x2rad(o,err);
            radstd = circ_std(merr,[],[],3);
            std = rad2x(o,radstd);
            this = strcat(o.name, 'Std Dev');
            imagesc(std,clim);
            title(this);
        end
    end
    
    methods (Static)
        function prms = defaultVMprms(type)
            switch upper(type)
                case 'NET2NET'
                    prms.k = 1;
                    prms.sigma = 0.37;
                    prms.v = 0;
                    prms.dim = 1;
                    
                case 'WORLD2NET'
                    prms.k = 20;
                    prms.sigma = 0.40;
                    prms.v = 0;
                    prms.dim = 1;
            end  
        end
    end
end