% CBSfunc.m
% Arthur Lee. Modified Jan.20.2020
% Function for calculating AUC of CBS and yhat from x using numerical approximation.

function out = CBSfunc(xpos,ypos,x,javaswitch)
% if javaswitch is true or missing, we call java function, which should be about 8~9 times faster than the matlab implementation
if nargin<3
    out = CBSAUC(xpos,ypos); % only xpos and ypos provided. calculating AUC
else
    if nargin < 4
        javaswitch = true;
    end
    x = x(:); % make x a column vector
    if javaswitch
        try
            out = javaMethod('getyhat','CBScalc',xpos,ypos,x); % try to call java function
        catch
            % adding java classpath
            javaaddpath(fileparts(which('CBSfunc'))); % assume that the java class file is in the same location as this file
            out = javaMethod('getyhat','CBScalc',xpos,ypos,x); % try to call java function
        end
    else
        out = CBScalc_slow(xpos,ypos,x);
    end
end
end

function AUC = CBSAUC(xpos,ypos)
AUC = 0;
partialAUC = @(x1,x2,x3,x4,y1,y2,y3,y4) (6*x2*y1-6*x1*y2-10*x1*y1-3*x1*y3+3*x3*y1-x1*y4-3*x2*y3+3*x3*y2+x4*y1-3*x2*y4+3*x4*y2-6*x3*y4+6*x4*y3+10*x4*y4)/20;
for i = 1:3:(length(xpos)-1)
    AUC = AUC + partialAUC(xpos(i),xpos(i+1),xpos(i+2),xpos(i+3),ypos(i),ypos(i+1),ypos(i+2),ypos(i+3));
end
end

function y = CBScalc_slow(xpos,ypos,x)
y = nan(size(x));
for i = 4:3:length(xpos)
    if i == 4
        idx = x<=xpos(i);
    else
        idx = xpos(i-3) < x & x <=xpos(i);
    end
    % calculate cubic equation coefficients
    a = -xpos(i-3)+3*xpos(i-2)-3*xpos(i-1)+xpos(i);
    b = 3*xpos(i-3)-6*xpos(i-2)+3*xpos(i-1);
    c = -3*xpos(i-3)+3*xpos(i-2);
    d = xpos(i-3)-x(idx);
    
    t = 0.5; % initial point
    ft = a*t^3 + b*t^2 + c*t+d; % f(t)
    delta = 0.25;
    
    while(max(abs(ft))>0.0000001) % using bisection, because other methods kept finding root outside of 0 1 (e.g., Newton-Raphson, Halley's Method)
        t = t-sign(ft).*delta;
        delta = delta/2;
        ft = a.*t.^3 + b.*t.^2 + c.*t+d;
    end
    
    % calculate using deCasteljau's algorithm for numerical stability
    y(idx) = (1-t).*((1-t).*((1-t).*ypos(i-3)+t.*ypos(i-2))+t.*((1-t).*ypos(i-2)+t.*ypos(i-1)))+t.*((1-t).*((1-t).*ypos(i-2)+t.*ypos(i-1))+t.*((1-t).*ypos(i-1)+t.*ypos(i)));
end
end
