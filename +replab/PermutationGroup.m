classdef PermutationGroup < replab.bsgs.BSGSGroup
% Represents a permutation group

    properties (SetAccess = protected)
        domainSize; % d when this group acts on {1, ..., d}
    end
    
    methods
                
        function self = PermutationGroup(domainSize, generators, orderOpt)
            if nargin < 3
                orderOpt = [];
            end
            A = replab.Permutations(domainSize).naturalAction;
            self@replab.bsgs.BSGSGroup(A, generators, orderOpt)
            self.domainSize = domainSize;
        end
        
        function R = naturalRepresentation(self)
            field = 'R15'; % TODO generalize
            hom = @(g) replab.Permutations(self.domainSize).toMatrix(g);
            images = cell(1, self.nGenerators);
            for i = 1:self.nGenerators
                images{i} = hom(self.generators{i});
            end
            target = replab.GeneralLinearGroup(self.domainSize, field);
            R = replab.FiniteGroupRep(self, images, target);
        end
        
    end
    
    methods (Static) % Factory methods
        
        function G = trivial(domainSize)
        % Constructs the trivial group on the domain
            G = replab.PermutationGroup(domainSize, {}, vpi(1));
        end
        
        function G = symmetric(n)
        % Constructs the symmetric group acting on n elements
            order = vpi(n);
            for i = 2:n-1
                order = order * vpi(i);
            end
            if n == 1
                G = replab.PermutationGroup.trivial(n);
            elseif n == 2
                G = replab.PermutationGroup(2, {[2 1]}, order);
            else
                G = replab.PermutationGroup(n, {[2:n 1] [2 1 3:n]}, order);
            end
        end
        
        function G = cyclic(n)
        % Constructs the cyclic group acting on n elements
            if n == 1
                G = replab.PermutationGroup.trivial(n);
            else
                G = replab.PermutationGroup(n, {[2:n 1]}, vpi(n));
            end
        end
        
        function G = fromGeneratorsAsMatrix(genMat)
            errmsg = 'To create the trivial group, pass zeros(0, n) so that the domain size n is known';
            domainSize = size(genMat, 2);
            assert(domainSize > 0, errmsg);
            nG = size(genMat, 1);
            generators = cell(1, nG);
            for i = 1:nG
                generators{i} = genMat(i, :);
            end
            G = replab.PermutationGroup(domainSize, generators);
        end
        
        function G = fromDomainSizeAndGenerators(domainSize, varargin)
            generators = varargin;
            G = replab.PermutationGroup(domainSize, generators);
        end
        
        function G = fromGenerators(varargin)
            errmsg = 'Argument list cannot be empty. Use fromDomainSizeAndGenerators instead';
            assert(nargin > 0, errmsg);
            g = varargin{1};
            domainSize = length(g);
            G = replab.PermutationGroup.fromDomainSizeAndGenerators(domainSize, varargin{:});
        end
        
    end
    
end
