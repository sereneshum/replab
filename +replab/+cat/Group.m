classdef Group
% Defines a group
    
    properties (Abstract, SetAccess = private)
        identity;
    end
    
    methods
        function xInv = inverse(self, x)
            error('Not implemented');
        end
            function z = compose(self, x, y, varargin)
            error('Not implemented');
        end
        function b = isIdentity(self, x)
            error('Not implemented');
        end
    end
    
    methods
        
        function x = evaluateWord(self, word, generators)
            x = self.identity;
            for i = 1:length(word.indices)
                g = generators{word.indices(i)};
                e = word.exponents(i);
                we = self.composeN(g, e);
                x = self.compose(x, we);
            end
        end
        
        function x = conjugate(self, by, on)
        % Returns "on" conjugated by "by", i.e.
        % by * on * by.inverse
            x = self.composeWithInverse(self.compose(by, on), by);
        end

        function z = composeWithInverse(self, x, y)
            z = self.compose(x, self.inverse(y));
        end
        
        function z = composeMany(self, varargin)
            if length(varargin) == 0
                z = self.identity;
            else
                z = varargin{1};
                for i = 2:length(varargin)
                    z = self.compose(z, varargin{i});
                end
            end
        end
        
        function y = composeN(self, x, n)
        % Computes y = x^e by repeated squaring
            if n < 0
                y = self.composeN(self.inverse(x), -n);
            elseif n == 0
                y = self.identity;
            else
                y = self.identity;
                while n > 1
                    if mod(n, 2) == 0 % n even
                        n = n / 2;
                    else % n odd
                        y = self.compose(x, y);
                        n = (n - 1)/2;
                    end
                    x = self.compose(x, x);
                end
                y = self.compose(x, y);
            end
        end
        
    end
    
end
