function  ws = windspeed(t, t1,t2, windtype, varargin)
%-------------------------------------------------------------------
%  generate a sin-wave wind gust at certain time range [t1, t2]
%  
%-------------------------------------------------------------------
if nargin > 4
    bias = varargin{1};
else
    bias = 0;
end

if nargin > 5
    amp = varargin{2};
else
    amp = 1;
end

switch windtype
    case 'step';
        ws = bias*(0<=t & t<t1)+...
            (amp+bias)*(t>=t1 & t<=t2)+...
            bias*(t>t2);
    case 'sine';
        ws = bias*(0<=t & t<t1)+...
            (amp*sin(pi/(t2-t1)*(t-t1))+bias)*(t>=t1 & t<=t2)+...
            bias*(t>t2);
    otherwise;
        ws = bias*(0<=t & t<t1)+...
            (amp*sin(pi/(t2-t1)*(t-t1))+bias)*(t>=t1 & t<=t2)...
            +bias*(t>t2);
end

end