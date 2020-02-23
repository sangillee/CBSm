% CBS_RC.m
%
% Fit either a 1-piece or 2-piece CBS latent utility function to binary risky choice data.
%
% The input data has n choices (ideally n > 100) between two reward options.
% Option 1 is receiving 'Amt1' with probability 'Prob1' and Option 2 is receiving 'Amt2' with probability 'Prob2' (e.g., $40 with 53% chance vs. $20 with 90% chance).
% One of the two options may be certain (i.e., prob = 1; e.g., $40 with 53% chance vs. $20 for sure).
% 'choice' should be 1 if option 1 is chosen, 0 if option 2 is chosen.
%
% 'choice Vector of 0s and 1s. 1 if the choice was option 1, 0 if the choice was option 2.
% 'Amt1' : Vector of positive real numbers. Reward amount of choice 1.
% 'Prob1' : Vector of positive real numbers between 0 and 1. Probability of winning the reward of choice 1.
% 'Amt2' : Vector of positive real numbers. Reward amount of choice 2.
% 'Prob2' : Vector of positive real numbers between 0 and 1. Probability of winning the reward of choice 2.
% 'numpiece' : Either 1 or 2. Number of CBS pieces to use.
% 'out' : A struct containing the following:
%       'type' : either 'CBS1' or 'CBS2' depending on the number of pieces
%       'LL' : log likelihood of the model
%       'numparam' : number of total parameters in the model
%       'scale' : scaling factor of the logit model
%       'xpos' : x coordinates of the fitted CBS function
%       'ypos' : y coordinates of the fitted CBS function
%       'AUC' : area under the curve of the fitted CBS function. Normalized to be between 0 and 1.

function out = CBS_RC(choice,Amt1,Prob1,Amt2,Prob2,numpiece)
% error checking
assert(all(choice == 0 | choice == 1),'Choice should be a vector of 0 or 1')
assert(all(Amt1 >= 0 & Amt2 >=0),'Negative amounts are not allowed')
assert(all(Prob1 >= 0 & Prob2 >= 0),'Negative delays are not allowed')
assert(numpiece == 1 || numpiece ==2,'Sorry! Only 1-piece and 2-piece CBS functions are supported at the moment.')

minpad = 1e-04; maxpad = 1-minpad; % pad around bounds because the solving algorithm tests values around the bounds
numfit = 40*numpiece; % number of search starting points

% checking to make sure probability is within 0 and 1
assert(all(Prob1<=1 | Prob2<=1),'prob not within [0 1]');

% parameter bounds and constraints
lb = [-36,minpad.*ones(1,6*numpiece-2)]; ub = [36,maxpad.*ones(1,6*numpiece-2)];
if numpiece==1 % active parameters (5): logbeta, x2, x3, y2, y3
    A = []; b = []; % no linear constraints
    nonlcon = []; % no non-linear constraints
    options = optimset('Display','off','Algorithm','sqp','TolCon',minpad);
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
    b = -minpad.*ones(8,1);
    nonlcon = @twopiece_nonlincon; % non-linear constraints
    options = optimset('Display','off','GradConstr','on','Algorithm','sqp','TolCon',minpad);
    x0 = [0,1/6,2/6,3/6,4/6,5/6, 1/6,2/6,3/6,4/6,5/6];
end

% fitting
problem = createOptimProblem('fmincon','x0',x0,'objective',@(x)negLL(x,Amt1,Prob1,Amt2,Prob2,choice),'lb',lb,'ub',ub,'Aineq',A,'bineq',b,'nonlcon',nonlcon,'options',options);
[outparam,fval] = run(MultiStart,problem,numfit);

% organizing output
out.type = ['CBS',num2str(numpiece)]; out.LL = -fval*length(choice);
out.numparam = size(outparam,2); out.scale = exp(outparam(1));
out.xpos = [0,outparam(2:((out.numparam+1)/2)),1]; out.ypos = [0,outparam(((out.numparam+1)/2+1):end),1];
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