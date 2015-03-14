classdef AnimalRegistry < handle
   properties (SetAccess = 'private')
       Animals = Animal.empty;
   end
   
   methods
       function addAnimal(obj,animal)
           if ~any(isequal(obj.Animals,animal))
               obj.Animals = [obj.Animals animal];
           end
       end
       
   end
end