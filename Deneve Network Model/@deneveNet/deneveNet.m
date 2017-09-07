classdef deneveNet < dynamicprops
    %A Deneve network model, in which one or more interconnected layers of neurons
    %represent population codes for sensory, motor, or cognitive variables
    %(e.g. position, motion, orientation, etc.). Layer-to-layer connections
    %can be feedforward, feedback, and any combination thereof across
    %the network.
    %
    %See deneveDemo.m for an example of how to define and connect network layers
    %and to run simulations.
    %
    %The user can gain control at different times during a simulation by
    %assigning callback functions to events (see constructor). Each function
    %should accept a deneveNet object as the sole argument.
    %
    %The state of the network can be logged automatically. This is done
    %within each layer object through the deneveLayer class.
    %
    % See Deneve, Latham, & Pouget (2001), "Efficient computation and cuintegration with noisy population codes."
    properties
        plotIt@logical = false;
        evtFun@struct = struct('preSim',@preSim,'timeZero',@timeZero,'beforeUpdate',@beforeUpdate, 'afterUpdate',@afterUpdate,'plot',@plot); %Structure of function handles for each event
    end
    
    properties (SetAccess = protected)
        name@char;
        layers@deneveLayer;
        
        %Simulation properties
       t                           %Time index for the current simulation
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
            
            %Assign the default event handlers
            allEvts = fieldnames(o.evtFun);
            for i=1:numel(allEvts)
                o.evtFun.(allEvts{i}) = eval(['@',allEvts{i}]);
            end
            
        end
        
        function addLayer(o,layer)
            %Add a deneveLayer to the network
            addprop(o,layer.name);
            o.(layer.name) = layer;
            
            %Store the handle
            o.layers = horzcat(o.layers,layer);
            
            %Register that I am the network parent of the layer
            layer.n = o;
        end
        
        function o = reset(o)
            %Reset the network responses to zero between simulations
            for i=1:o.nLayers
                o.layers(i).reset();
            end
        end
        
        function o = noisy(o,state)
            %Convenience function to switch noise on (true) or off (false) for all layers of the network
            for i=1:o.nLayers
                o.layers(i).noisy = state;
            end
        end
        
        function o = update(o)
            %Update all of the layers
            
            %Call timunsonerep on all layers
%             for i=1:o.nLayers
%                 o.layers(i).bufferResp(); %this becomes t minus one
%             end
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
        
        function o = logState(o,newLog)
            %Log the state of all layers
            for i=1:o.nLayers
                o.layers(i).logState('newLog',newLog);
            end
        end
        
        function run(o,varargin)
            %Runs a simulation for the network in its current state.
            %All setting of input weights, initial states etc. should be
            %done outside this function. Done this way in the hope that it
            %will be compatible with a parallel approach (parfor).
            p = inputParser;
            p.addParameter('nIter',10,@(x) isnumeric(x) & x > 0);   %How many network updates should we run?
            p.addParameter('logState',true,@islogical);
            p.parse(varargin{:});
            p = p.Results;
            
            %Run pre-sim code (e.g. to set the initial state of some layers)
            o.evtFun.preSim(o);
            
            %Run the simulation
            for tInd=1:p.nIter
                
                %Set the time
                o.t = tInd;
                
                %Run time-zero function
                isTimeZero = tInd==1;
                if isTimeZero, o.evtFun.timeZero(o); end
                
                %Run pre-update code
                o.evtFun.beforeUpdate(o);
                
                %Calculate the new state of the network
                o.update();
                
                %Run post-update code
                o.evtFun.afterUpdate(o);
                
                %Run plotting code
                if o.plotIt, o.evtFun.plot(o); end
                
                %Log the current state of the network
                if p.logState
                    o.logState(isTimeZero);
                end
            end
        end
        
        function preSim(o)
            
        end
        
        function timeZero(o)
            %To be overloaded in child sims
        end
        
        function beforeUpdate(o)
            
        end
        
        function afterUpdate(o)
            
        end
        
        function plot(o)
            
        end
    end
end