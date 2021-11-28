classdef PSDOptions < matlab.mixin.SetGet

    
    properties        
        nPointfilter =1 % factor that refines the accuracy calculating reduced dynamics on SSM
    end
    
    methods      
        function set.nPointfilter(obj,n)
            obj.nPointfilter = n;
        end  
        
    end
end

