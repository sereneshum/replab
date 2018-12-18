classdef FinitelyGeneratedGroup < handle
   
    properties
        generators; % 1 x nG cell array of generators
    end
    
    methods
                
        function n = nGenerators(self) 
        % Returns the number of generators
            n = length(self.generators);
        end

        function p = generator(self, i) 
        % Returns the i-th generator
            p = self.generators(i);
        end
                
        function p = generatorInverse(self, i)
        % Returns the inverse of the i-th generator of this group
            p = self.inverse(self.generator(i));
        end
        
        function b = isTrivial(self)
            b = self.nGenerators == 0;
        end
        
        function x = evaluateWord(self, word)
            x = self.G.identity;
            for i = 1:length(word.indices)
                g = self.generators{word.indices(i)};
                e = word.exponents(i);
                we = self.G.composeN(g, e);
                x = self.G.compose(x, we);
            end
        end
        
        function f = freeGroup(self)
            f = replab.FreeGroup(self.nGenerators);
        end
            
        function phi = morphismFromFreeGroup(self)
            phi = replab.cat.GroupMorphismFun(@(w) self.evaluateWord(word), self.freeGroup, self);
        end

    end
    
end
