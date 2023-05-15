function windowtitle(varargin)
% windowtitle(titlestr)
% windowtitle(fig, titlestr)
%
% Sets the title of the plot window to titlestr.
narginchk(1, 2);
if nargin() == 1
    titlestr = varargin{1};
    fig = gcf();
else
    titlestr = varargin{2};
    fig = varargin{1};
end
set(fig, 'NumberTitle', 'off', 'Name', titlestr);
end%function

