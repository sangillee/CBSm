% CBSfunc.m
%
% Calculate either the Area Under the Curve (AUC) of a CBS function, or calculate the y coordinates of CBS function given x.
% 'xpos' : Vector of real numbers of length 1+3n (n = 1, 2, 3, ...), corresponding to Bezier points' x-coordinates of a CBS function
% 'ypos' : Vector of real numbers of length 1+3n (n = 1, 2, 3, ...), corresponding to Bezier points' y-coordinates of a CBS function
% 'x': Vector of real numbers, corresponding to x-coordinates of a CBS function. Default value is Null.
% 'out' If x is provided, return y coordinates corresponding to x. If x is not provided, return AUC.
% 
% examples
% CBSfunc([0,0.3,0.6,1],[0.5, 0.2, 0.7, 0.9])
% CBSfunc([0,0.3,0.6,1],[0.5, 0.2, 0.7, 0.9],0:.01:1)

function out = CBSfunc(xpos,ypos,x)
assert(length(xpos)==length(ypos),'length of xpos and ypos different!')
assert(length(xpos) >= 4,'length of xpos and ypos too short. They must have at least 4 elements')
assert(mod(length(xpos),3) == 1,'unexpected length of xpos and ypos. They should be 3n+1 (n = 1, 2, ...)')
if nargin<3
    out = CBSAUC(xpos,ypos); % only xpos and ypos provided. calculating AUC
else
    x = x(:); % make x a column vector
    try
        out = javaMethod('getyhat','CBScalc',xpos,ypos,x); % try to call java function
    catch
        % adding java classpath
        javaaddpath(fileparts(which('CBSfunc'))); % assume that the java class file is in the same location as this file
        out = javaMethod('getyhat','CBScalc',xpos,ypos,x); % try to call java function
    end
end
end

function AUC = CBSAUC(xpos,ypos)
AUC = 0;
for i = 1:3:(length(xpos)-1)
    AUC = AUC + partialAUC(xpos(i),xpos(i+1),xpos(i+2),xpos(i+3),ypos(i),ypos(i+1),ypos(i+2),ypos(i+3));
end
end

function out = partialAUC(x1,x2,x3,x4,y1,y2,y3,y4)
check1 =  -sqrt((x4-x3)*(x2-x1)) < (x3-x2);
check2 = x1 <= x2;
check3 = x3 <= x4;
if ~check1 || ~check2 || ~check3
    print('CBS x coordinates not a monotonic function of t. Multiple y for x may exist. AUC may be inaccurate')
end
out = (6*x2*y1-6*x1*y2-10*x1*y1-3*x1*y3+3*x3*y1-x1*y4-3*x2*y3+3*x3*y2+x4*y1-3*x2*y4+3*x4*y2-6*x3*y4+6*x4*y3+10*x4*y4)/20;
end