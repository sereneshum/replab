function u = uminus(a)
% -  Unary minus.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

u = overload(mfilename, a);

% $Id: uminus.m 1004 2017-11-15 17:14:09Z sangwine $
