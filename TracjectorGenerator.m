
function yref = TracjectorGenerator(x,tType)
%-------------------------------------------------------------------------
%       1  -->  NLTV System -2, Trajectroy - Curve;
%       2 -->   NLTV System -1, Trajectroy - Curve;
%       3 -->   NLTV System -1, Trajectroy - Step;
%%
switch tType
    case 1;
        y1 = 5*sin(x/50) + 2*cos(x/20);
        y2 = 2*sin(x/50) + 5*cos(x/20);
    case 2;
       y1 = 0.75*sin(x*pi/8) + 0.5*cos(pi*x/4);
       y2 = 0.5*cos(x*pi/8) + 0.5*sin(pi*x/4);
    case 3;
        y1=0.4*(x>=0 & x<=500)+0.7*(x>500 & x<=1000)+0.5*(x>1000 & x<=1500);
        y2=0.6*(x>=0 & x<=300)+0.8*(x>300 & x<=700)...
            +0.7*(x>700 & x<=1200)+0.5*(x>1200 & x<=1500);
    case 4;
        y1=0.4*(x>=0 & x<=100)+0.7*(x>100 & x<=150)+0.5*(x>150 & x<=200);
        y2=0.6*(x>=0 & x<=80)+0.8*(x>80 & x<=130)...
            +0.7*(x>130 & x<=180)+0.5*(x>180 & x<=200);
end
yref = [y1', y2'];
return