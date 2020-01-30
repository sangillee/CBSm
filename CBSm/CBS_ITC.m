% CBS_ITC.m
%
% Fit either a 1-piece or 2-piece CBS latent utility function to binary intertemporal choice data.
%
% The input data has n choices (ideally n > 100) between two reward options.
% Option 1 is receiving 'Amt1' in 'Delay1' and Option 2 is receiving 'Amt2' in 'Delay2' (e.g., $40 in 20 days vs. $20 in 3 days).
% One of the two options may be immediate (i.e., delay = 0; e.g., $40 in 20 days vs. $20 today).
% 'choice' should be 1 if option 1 is chosen, 0 if option 2 is chosen.
%
% 'choice' : Vector of 0s and 1s. 1 if the choice was option 1, 0 if the choice was option 2.
% 'Amt1' : Vector of positive real numbers. Reward amount of choice 1.
% 'Delay1' : Vector of positive real numbers. Delay until the reward of choice 1.
% 'Amt2' : Vector of positive real numbers. Reward amount of choice 2.
% 'Delay2' : Vector of positive real numbers. Delay until the reward of choice 2.
% 'numpiece' : Either 1 or 2. Number of CBS pieces to use.
% 'out' : A struct containing the following:
%       'type' : either 'CBS1' or 'CBS2' depending on the number of pieces
%       'LL' : log likelihood of the model
%       'numparam' : number of total parameters in the model
%       'scale' : scaling factor of the logit model
%       'xpos': x coordinates of the fitted CBS function
%       'ypos': y coordinates of the fitted CBS function
%       'AUC': area under the curve of the fitted CBS function. Normalized to be between 0 and 1.
%       'normD' : The domain of CBS function runs from 0 to \code{normD}. Specifically, this is the constant used to normalize all delays between 0 and 1, since CBS is fitted in a unit square first and then scaled up.

function out = CBS_ITC(choice,Amt1,Delay1,Amt2,Delay2,numpiece)
% error checking
assert(all(choice == 0 | choice == 1),"Choice should be a vector of 0 or 1")
assert(all(Amt1 >= 0 & Amt2 >=0),"Negative amounts are not allowed")
assert(all(Delay1 >= 0 & Delay2 >= 0),"Negative delays are not allowed")
assert(numpiece == 1 || numpiece ==2,"Sorry! Only 1-piece and 2-piece CBS functions are supported at the moment.")

minpad = 1e-04; maxpad = 1-minpad; % pad around bounds because the solving algorithm tests values around the bounds
numfit = 40*numpiece; % number of search starting points

% normalizing delay to [0 1] for easier parameter search
nD = 1;
if (any(Delay1 > 1) || any(Delay2 > 1))
    nD = max([Delay1;Delay2]); Delay1 = Delay1./nD; Delay2 = Delay2./nD;
end

% parameter bounds and constraints
lb = [-36,minpad.*ones(1,6*numpiece-1)]; ub = [36,maxpad.*ones(1,6*numpiece-1)];
if numpiece==1 % active parameters (6): logbeta, x2, x3, y2, y3, y4
    A = [0,0,0,0,-1,1;0,0,0,-1,0,1]; b = [-minpad,-minpad]; % linear constraints: y4-y3<0, y4-y2<0
    nonlcon = []; % no non-linear constraints
    options = optimset('Display','off','Algorithm','sqp','TolCon',minpad);
    x0 = [0,1/3,2/3,2/3,1/3,0.01];
elseif numpiece == 2 % active parameters (12): logbeta, x2,x3,x4,x5,x6, y2,y3,y4,y5,y6,y7
    % linear constraints:
    A = [0,1,0,-1,0,0, zeros(1,6);... % x2-x4<0
        0,0,1,-1,0,0, zeros(1,6);...  % x3-x4<0
        0,0,0,1,-1,0, zeros(1,6);...  % x4-x5<0
        0,0,0,1,0,-1, zeros(1,6);...  % x4-x6<0
        zeros(1,6), -1,0,1,0,0,0;...  % y4-y2<0
        zeros(1,6), 0,-1,1,0,0,0;...  % y4-y3<0
        zeros(1,6), 0,0,-1,1,0,0;...  % y5-y4<0
        zeros(1,6), 0,0,-1,0,1,0;...  % y6-y4<0
        zeros(1,6), 0,0,0,-1,0,1;...  % y7-y5<0
        zeros(1,6), 0,0,0,0,-1,1];    % y7-y6<0
    b = -minpad.*ones(10,1);
    nonlcon = @twopiece_nonlincon; % non-linear constraints
    options = optimset('Display','off','GradConstr','on','Algorithm','sqp','TolCon',minpad);
    x0 = [0,1/6,2/6,3/6,4/6,5/6, 5/6,4/6,3/6,2/6,1/6,0.01];
end

% fitting
problem = createOptimProblem('fmincon','x0',x0,'objective',@(x)negLL(x,Amt1,Delay1,Amt2,Delay2,choice),'lb',lb,'ub',ub,'Aineq',A,'bineq',b,'nonlcon',nonlcon,'options',options);
[outparam,fval] = run(MultiStart,problem,numfit);

% organizing output
out.type = ['CBS',num2str(numpiece)]; out.LL = -fval*length(choice);
out.numparam = size(outparam,2); out.scale = exp(outparam(1));
out.xpos = nD.*[0,outparam(2:(out.numparam/2)),1]; out.ypos = [1,outparam((out.numparam/2 +1):end)];
out.AUC = CBSfunc(out.xpos,out.ypos); out.normD = nD;
end

function negLL = negLL(params,Amt1,Var1,Amt2,Var2,choice)
cutoff = size(params,2)/2;
yhat1 = CBSfunc([0,params(2:cutoff),1],[1,params((cutoff+1):end)],Var1);
yhat2 = CBSfunc([0,params(2:cutoff),1],[1,params((cutoff+1):end)],Var2);
DV = Amt1.*yhat1 - Amt2.*yhat2; % diff between utilities
DV(choice==0) = -DV(choice==0); % utility toward choice
reg = -exp(params(1)).*DV; % scaling by noise parameter
logp = -log(1+exp(reg)); % directly calculating logp
logp(reg>709) = -reg(reg>709); % log(realmax) is about 709.7827, reg>709 would have lost precision
negLL = -mean(logp); % making normalized LL space. I know this shouldn't be different from sum(logp), but this gives much reliable convergences...(for some weird reason)
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
gradc = [0,0,-2*x3+2*x4,-2*x4+2*x3,0,0,0,-2*y3+2*y4,-2*y4+2*y3,0,0,0;...
    0,0,0,-2*x4+2*x5,-2*x5+2*x4,0,0,0,-2*y4+2*y5,-2*y5+2*y4,0,0]';
gradceq = [0,0,y4-y5,y5-y3,y3-y4,0,0,x5-x4,x3-x5,x4-x3,0,0]';
end