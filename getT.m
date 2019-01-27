function [t]=getT(s,fs,varargin)
% Simple function that creates a time vector "t" with the same length as
% inputed vector "s" with sampling frequency "fs"
%


% Process parameters  %
plotflag='';

n = 0;
while n < length(varargin)
n = n + 1;
	if strcmp(varargin{n}, 'plot')
        plotflag = 'plot';
%	elseif strcmp(varargin{n}, ''
	end
end

%start code

t=linspace(0,(length(s)-1)/fs,length(s))';


if strcmp(plotflag,'plot')




end

end