classdef Animal < handle
    
    properties
        %CSCFiles
        ID
        Group = Group.empty;
        Sessions = Session.empty;
    end
    
    %% Get and Set Methods.  
    % See:  web([docroot '/techdoc/matlab_oop/brgsek9-1.html#brgsek9-3'])
    methods
        
        function set.Group(obj,group)
            % Error checking
            validateattributes(group,{'Group'},{});
            obj.Group = group;
            obj.Group.addAnimal(obj);
        end
    end
    methods
        function obj = Animal(varargin)
            % Animal Constructor:
                %       Animal(allAnimals,group);
                %          allAnimals is the AnimalsRegistry object that
                %          contains an array of all known animals
                %
                %          group is optional
            if nargin
                animalsarray = varargin{1};
                animalsarray.addAnimal(obj);
            end
            if nargin == 2
                group = varargin{2};
                obj.Group = group;
            end
        end
        function id = getGroupID(obj)
            id = obj.Group.ID;
        end
        function addSession(obj,session)
            obj.Sessions = [obj.Sessions session];
        end
        function trials = getTrials(obj,param,val)
            trials = [];
            for k = 1:length(obj)
                trials = [trials obj(k).Sessions.getTrials(param,val)]; %#ok<AGROW>
            end
        end
    end
    
    methods (Static)
        % Doc on static methods: web([docroot '/techdoc/matlab_oop/brdqiu3.html#brdsi37'])
        function animal = getAnimalFromID(animalsregistry,id)
            animal = animalsregistry.Animals(strcmp({animalsregistry.Animals.ID},id));
        end
    end
end
