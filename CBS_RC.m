function out = CBS_RC(choice,Amt1,Var1,Amt2,Var2,numpiece)
% function for using fmincon to fit CBS 1-piece or 2-piece functions to risky choice
% choice would be 1 if they chose option 1 or 0 if they chose option 2
% Var should be probability
assert( all(Var1<=1) && all(Var2<=1) && all(Var1>=0) && all(Var2>=0) ,'prob not within [0 1]');

lb = [-36,zeros(1,6*numpiece-2)]; % lower bounds
ub = [36,ones(1,6*numpiece-2)]; % upper bounds
numfit = 40*numpiece; % number of search starting points

if numpiece==1 % active parameters (5): logbeta, x2, x3, y2, y3
    A = []; b = []; % no linear constraints
    nonlcon = []; % no non-linear constraints
    options = optimset('Display','off','Algorithm','sqp');
    x0 = [0,1/3,2/3,1/3,2/3];
elseif numpiece == 2 % active parameters (11): logbeta, x2,x3,x4,x5,x6, y2,y3,y4,y5,y6
    % linear constraints:
    A = [0,1,0,-1,0,0, zeros(1,5);...   % x2-x4<0
         0,0,1,-1,0,0, zeros(1,5);...   % x3-x4<0
         0,0,0,1,-1,0, zeros(1,5);...   % x4-x5<0
         0,0,0,1,0,-1, zeros(1,5);...   % x4-x6<0
         zeros(1,6), 1,0,-1,0,0;...     % y2-y4<0
         zeros(1,6), 0,1,-1,0,0;...     % y3-y4<0
         zeros(1,6), 0,0,1,-1,0;...     % y4-y5<0
         zeros(1,6), 0,0,1,0,-1];       % y4-y6<0
    
    
    b = zeros(8,1);
    nonlcon = @twopiece_nonlincon; % non-linear constraints
    options = optimset('Display','off','GradConstr','on','Algorithm','sqp');
    x0 = [0,1/6,2/6,3/6,4/6,5/6, 1/6,2/6,3/6,4/6,5/6];
end
problem = createOptimProblem('fmincon','x0',x0,'objective',@(x)negLL(x,Amt1,Var1,Amt2,Var2,choice),'lb',lb,'ub',ub,'Aineq',A,'bineq',b,'nonlcon',nonlcon,'options',options);
[outparam,fval] = run(MultiStart,problem,numfit);
out.type = ['CBS',num2str(numpiece)];
out.LL = -fval*length(choice);
out.numparam = size(outparam,2);
out.scale = exp(outparam(1));
out.xpos = [0,outparam(2:((out.numparam+1)/2)),1];
out.ypos = [0,outparam(((out.numparam+1)/2+1):end),1];
out.AUC = CBSfunc(out.xpos,out.ypos);
end

function negLL = negLL(params,Amt1,Var1,Amt2,Var2,choice)
cutoff = (size(params,2)+1)/2;
yhat1 = CBSfunc([0,params(2:cutoff),1],[0,params((cutoff+1):end),1],Var1);
yhat2 = CBSfunc([0,params(2:cutoff),1],[0,params((cutoff+1):end),1],Var2);
DV = Amt1.*yhat1 - Amt2.*yhat2; % diff between utilities
DV(choice==0) = -DV(choice==0); % utility toward choice
reg = -exp(params(1)).*DV; % scaling by noise parameter
logp = -log(1+exp(reg)); % directly calculating logp
logp(reg>709) = -reg(reg>709); % log(realmax) is about 709.7827, reg>709 would have lost precision
negLL = -mean(logp); % making normalized LL
end

function [c,ceq,gradc,gradceq] = twopiece_nonlincon(params)
minhandle = 0.1;
x3 = params(3); x4 = params(4); x5 = params(5); y3 = params(8); y4 = params(9); y5 = params(10);
% non-linear inequalities: 0.1^2 -(x4-x3)^2 -(y4-y3)^2 < 0, 0.1^2-(x5-x4)^2-(y5-y4)^2 < 0
c(1) = minhandle^2 -(x4-x3)^2 -(y4-y3)^2;
c(2) = minhandle^2 -(x5-x4)^2 -(y5-y4)^2;
% non-linear equalities: (x4-x3)/(y4-y3) = (x5-x4)/(y5-y4)
ceq = (x4-x3)*(y5-y4)-(x5-x4)*(y4-y3);
% calculating gradients for non-linear inequalities
gradc = [0,0,-2*x3+2*x4,-2*x4+2*x3,0,0,0,-2*y3+2*y4,-2*y4+2*y3,0,0;...
    0,0,0,-2*x4+2*x5,-2*x5+2*x4,0,0,0,-2*y4+2*y5,-2*y5+2*y4,0]';
gradceq = [0,0,y4-y5,y5-y3,y3-y4,0,0,x5-x4,x3-x5,x4-x3,0]';
end