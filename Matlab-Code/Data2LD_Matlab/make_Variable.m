function variableStr = make_Variable(name, order, XCell, FCell)
%  make_Variable assembles four arguments into a struct object that is
%  used by function make_Model  to set up a linear dynamicsystem object.
%  Arguments are:
%  NAME  ... A string to be used as the name of the variable
%  ORDER ... The order of the derivative on the left side
%  XCELL ... A cell array containing specifications for the 
%            homogeneous terms.  See make_Xterm for more details.
%  FCELL ... A cell array containing specifications for tthe 
%            forcing terms.  See make_Fterm for more details.

%  Last modified 26 January 2019

if nargin < 4, FCell  = [];  end
if nargin < 3, XCell  = [];  end
if nargin < 2, order  = 1;   end
if nargin < 1, name   = '';  end
if floor(order) ~= order || order < 1
    error('Argument ORDER is not a positive integer.');
end
variableStr.name  = name;
variableStr.order = order;
variableStr.XCell = XCell;
variableStr.FCell = FCell;
end
