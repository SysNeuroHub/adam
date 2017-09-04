classdef deneveSim
    %Class to run a simulation using a deneveNet object.
    %Scripts can gain control at different times during the
    %simulation by assigning a callback function to different
    %events (see constructor). Each function should accept a sim object (s) as
    %the first argument).
    properties
        evtFun@struct = struct('preSim',@preSim,'timeZero',@timeZero,'beforeUpdate',@beforeUpdate, 'afterUpdate',@afterUpdate,'plot',@plot); %Structure of function handles for each event
        t                           %Time index for the current simulation
        plotIt@logical = false;
%         nSims
%         sInd                        %Simulation index
    end
    
    properties (SetAccess = protected)
        net@deneveNet;              %The deneveNet object used for this simulation
        log;
    end
    
    methods
        function s = deneveSim(n)
            %Constructor.
            s.net = n;
            
            %Assign the default event handlers
            allEvts = fieldnames(s.evtFun);
            for i=1:numel(allEvts)
                s.evtFun.(allEvts{i}) = eval(['@',allEvts{i}]);
            end
        end
        
        function run(s,varargin)
            %Runs a simulation for the network in its current state.
            %All setting of input weights, initial states etc. should be
            %done outside this function. Done this way in the hope that it
            %will be compatible with a parallel approach (parfor).
            p = inputParser;
            p.addParameter('nIter',10,@(x) isnumeric(x) & x > 0); %How many network updates should we run?
            p.parse(varargin{:});
            p = p.Results;

            %Allocate memory to log network state
            allocLog(s.net,s.nIter,1); %Allocates for all layers. Call allocLog on deneveLayer objects directly if you don't want all layers to log
            
            %Run pre-sim code (e.g. to set the initial state of some layers)
            s.evtFun.preSim(s);
            
            %Run the simulation
            for tInd=1:p.nIter
                                           
                %Set the time
                s.t = tInd;
                
                %Run time-zero function
                if tInd==1, s.evtFun.timeZero(s); end
 
                %Run pre-update code
                s.evtFun.beforeUpdate(s);
                
                %Calculate the new state of the network
                s.net.update();

                %Run post-update code
                s.evtFun.afterUpdate(s);
                
                %Run plotting code
                if s.plotIt, s.evtFun.plot(s); end

                %Probably should be some logging here (e.g. stim used for this sim at this time?)
            end
        end
        
        function preSim(s)
            
        end
        
        function timeZero(s)
            %To be overloaded in child sims
        end
 
        function beforeUpdate(s)
            
        end
        
        function afterUpdate(s)
            
        end
        
        function plot(s)

        end   
    end
end