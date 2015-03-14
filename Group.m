classdef Group < handle
    properties
        Animals
        ID = 5;
    end
    
    methods
        function addAnimal(obj,animal)
            if ~any(obj.Animals==animal)
                obj.Animals = [obj.Animals animal];
            end
        end
        function trials = getTrials(obj,param,val)
            trials = obj.Animals.getTrials(param,val);
        end
    end
end