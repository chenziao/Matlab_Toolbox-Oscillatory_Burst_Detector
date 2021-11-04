%
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YTEA101
% Project Title: Particle Swarm Optimization Video Tutorial
% Publisher: Yarpiz (www.yarpiz.com)
%
% Developer and Instructor: S. Mostapha Kalami Heris (Member of Yarpiz Team)
%
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function out = PSO(problem, params)

%% Problem Definiton
CostFunction = problem.CostFunction;  % Cost Function
if isfield(problem,'LUbounds') && ~isempty(problem.LUbounds)
    if size(problem.LUbounds,2)~=2
        error('2 columns are required for "LUbounds".');
    end
    nVar = size(problem.LUbounds,1);    % Number of Decision Variables
    VarMin = problem.LUbounds(:,1)';    % Lower Bound of Decision Variables
    VarMax = problem.LUbounds(:,2)';    % Upper Bound of Decision Variables
else
    nVar = problem.nVar;        % Number of Unknown (Decision) Variables
    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables
    if isscalar(VarMin),    VarMin = repmat(VarMin,[1,nVar]);
    else	VarMin = reshape(VarMin,[1,nVar]);	end
    if isscalar(VarMax),    VarMax = repmat(VarMax,[1,nVar]);
    else	VarMax = reshape(VarMax,[1,nVar]);	end
end

%% Parameters of PSO
MaxIt = params.MaxIt;   % Maximum Number of Iterations
nPop = params.nPop;     % Population Size (Swarm Size)
w = params.w;           % Intertia Coefficient
wdamp = params.wdamp;   % Damping Ratio of Inertia Coefficient
c1 = params.c1;         % Personal Acceleration Coefficient
c2 = params.c2;         % Social Acceleration Coefficient
withNeighbor = params.withNeighbor;

% The Flag for Showing Iteration Information
ShowIterInfo = params.ShowIterInfo;
Visualize = params.Visualize;
if Visualize
    if isfield(params,'VarDisp') && numel(params.VarDisp)==2
        VarDisp = params.VarDisp;
    else
        if nVar>1,	VarDisp = [1,2];	else	VarDisp = [1,0];	end
    end
end

if isfield(params,'MinIt') && isscalar(params.MinIt)
    MinIt = params.MinIt;   % Minimum Number of Iterations
else
    MinIt = 0;
end

if isfield(params,'MaxStall') && isscalar(params.MaxStall)
    MaxStall = params.MaxStall;
else
    MaxStall = MaxIt;
end
if isfield(params,'FunTol') && isscalar(params.FunTol)
    FunTol = max(params.FunTol,0);
else
    FunTol = 0;
end

if withNeighbor
    if isfield(params,'NeighborPeriod') && isscalar(params.NeighborPeriod)...
            && params.NeighborPeriod>0 && params.NeighborPeriod<=MaxIt
        NbrPd = ceil(params.NeighborPeriod);
    else
        NbrPd = MaxIt;
    end
    if ~isfield(params,'NeighborSize') || isempty(params.NeighborSize)
        nNbr = (nPop-1)*[1,1];
    else
        nNbr = zeros(1,numel(params.NeighborSize));
        for i = 1:numel(params.NeighborSize)
            if params.NeighborSize(i)<=1
                nNbr(i) = max(round((nPop-1)*params.NeighborSize(i)),2);
            else
                nNbr(i) = max(min(round(params.NeighborSize(i)),nPop-1),2);
            end
        end
    end
    nupdate = ceil(MaxIt/NbrPd);
    updateNbr = reshape([round(linspace(nNbr(1),nNbr(end),nupdate));...
        repmat(zeros(NbrPd-1,1),[1,nupdate])],1,[]);
end

if isfield(params,'InitVel') && isscalar(params.InitVel) && params.InitVel
    InitVel = true;
else
    InitVel = false;
end
MaxVelocity = (VarMax-VarMin)/4;
MinVelocity = -MaxVelocity;

if isfield(params,'randseed') && isscalar(params.randseed) && params.randseed>0
    rng(params.randseed);
end

%% Initialization
% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
if withNeighbor,    empty_particle.LocalBest = empty_particle.Best; end
% Create Population Array
particle = repmat(empty_particle, nPop, 1);

% Initialize Population Members
% Generate Random Solution
randval = cellfun(@(x) (VarMin+VarMax)/2+(VarMax-VarMin).*...
    (rand(1,nVar)-0.5), cell(nPop,1),'UniformOutput',false);
[particle.Position] = deal(randval{:});
% Initialize Velocity
if InitVel
    randval = mat2cell(bsxfun(@times,2*rand(nPop,nVar)-1,MaxVelocity),ones(1,nPop));
    [particle.Velocity] = deal(randval{:});
else
    [particle.Velocity] = deal(zeros([1,nVar]));
end
randstate = rng;	% Store random state temporally
for i=1:nPop
    % Evaluation
    particle(i).Cost = CostFunction(particle(i).Position);
    % Update the Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
end
% Array to Hold Best Cost Value on Each Iteration
BestCosts = NaN(MaxIt+1,1);
% Initialize Global Best
[BestCosts(1),Ibest] = min([particle.Cost]);
GlobalBest = particle(Ibest).Best;
HistoryBest = BestCosts;

% Initialize Number of Stall Iterations
numStall = 0;

message = num2str([BestCosts(1),HistoryBest(1)],...
    'Iteration 0: Best Cost = %f, History Best = %f, Stall Iterations = 0');
if ShowIterInfo,	fprintf(message);	msglen = numel(message);	end
if Visualize,	fig = plotparticles;	end

%% Main Loop of PSO
for it = 1:MaxIt
    rng(randstate);     % Restore random state
    % Update Neighborhood and Local Best
    if withNeighbor
        if updateNbr(it),	Neighbors = genNeighbors(updateNbr(it));	end
        updateLocalBest(arrayfun(@(x) x.Best.Cost,particle));
    end
    % Generate Random Vectors for Velocity
    u = rand(2,nVar,nPop);
    randstate = rng;	% Store random state temporally
    for i=1:nPop
        % Update Velocity
        if withNeighbor,	SocialPos = particle(i).LocalBest.Position;
        else	SocialPos = GlobalBest.Position;	end
        particle(i).Velocity = w*particle(i).Velocity ...
            + c1*u(1,:,i).*(particle(i).Best.Position - particle(i).Position) ...
            + c2*u(2,:,i).*(SocialPos - particle(i).Position);
        % Apply Velocity Limits
        particle(i).Velocity = min(max(particle(i).Velocity,MinVelocity),MaxVelocity);
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        % Apply Lower and Upper Bound Limits
        particle(i).Position = min(max(particle(i).Position,VarMin),VarMax);
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        % Update Personal Best
        if particle(i).Cost < particle(i).Best.Cost
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
        end
    end
    
    % Update Global Best
    [BestCosts(it+1),Ibest] = min([particle.Cost]);
    % Store the Best Cost Value
    GlobalChange = HistoryBest(it)-BestCosts(it+1);
    if GlobalChange>0
        GlobalBest = particle(Ibest).Best;
        HistoryBest(it+1) = BestCosts(it+1);
    else
        HistoryBest(it+1) = HistoryBest(it);
    end
    % Update stall number
    if GlobalChange<=FunTol
        numStall = numStall+1;
    else
        numStall = 0;
    end
    
    % Damping Inertia Coefficient
    w = w * wdamp;
    
    % Display Iteration Information
    message = [num2str(it,'Iteration %d'),...
        num2str(BestCosts(it+1),': Best Cost = %f'),...
        num2str(HistoryBest(it+1),', History Best = %f'),...
        num2str(numStall,', Stall Iterations = %d')];
    if withNeighbor
        message = [message,num2str(size(Neighbors,2),', Neighbor Size = %d')];
    end
    if ShowIterInfo==1
        fprintf([repmat('\b',1,msglen),message]);
        msglen = numel(message);
    end
    if Visualize,	plotparticles(fig);	end
    
    % Check termination criteria
    if it>=MinIt && numStall>=MaxStall,	out = output(1);	return;	end
end

out = output(-1);

%% Nested Functions
    function Nbrs = genNeighbors(n)
        Nbrs = zeros(nPop,n);
        for j = 1:nPop
            iNbr = randperm(nPop-1,n);
            Nbrs(j,:) = [iNbr(iNbr<j),iNbr(iNbr>=j)+1];
        end
    end

    function updateLocalBest(SelfBestCosts)
        for j = 1:nPop
            [~,Ilbest] = min(SelfBestCosts(Neighbors(j,:)));
            particle(j).LocalBest = particle(Neighbors(j,Ilbest)).Best;
        end
    end

    function fig = plotparticles(fig)
        if nargin<1
            fig = figure;   grid on;	withNbr = false;    it=0;
        else
            figure(fig);
            if ~VarDisp(2),	ylm = get(gca,'ylim');	end
            clf;    withNbr = withNeighbor;
        end
        subplot(211);   hold on;
        x = arrayfun(@(x) x.Position(VarDisp(1)),particle);
        xs = arrayfun(@(x) x.Best.Position(VarDisp(1)),particle);
        xg = GlobalBest.Position(VarDisp(1));
        if VarDisp(2)
            y = arrayfun(@(x) x.Position(VarDisp(2)),particle);
            ys = arrayfun(@(x) x.Best.Position(VarDisp(2)),particle);
            yg = GlobalBest.Position(VarDisp(2));
        else
            y = [particle.Cost];    yg = GlobalBest.Cost;
            ys = arrayfun(@(x) x.Best.Cost,particle);
        end
        for j = 1:nPop
            plot([xs(j),x(j)],[ys(j),y(j)],'b.-');
        end
        plot(xs,ys,'m.','MarkerSize',7);
        if withNbr
            xl = arrayfun(@(x) x.LocalBest.Position(VarDisp(1)),particle);
            if VarDisp(2)
                yl= arrayfun(@(x) x.LocalBest.Position(VarDisp(2)),particle);
            else
                yl = arrayfun(@(x) x.LocalBest.Cost,particle);
            end
            plot(xl,yl,'co','MarkerSize',5)
        end
        plot(xg,yg,'ro','MarkerSize',8);
        xlim([VarMin(VarDisp(1)),VarMax(VarDisp(1))]);
        xlabel(num2str(VarDisp(1),'Var%d'));
        if VarDisp(2)
            ylim([VarMin(VarDisp(2)),VarMax(VarDisp(2))]);
            ylabel(num2str(VarDisp(2),'Var%d'));
        else
            if exist('yl','var')
                ylnew = get(gca,'ylim');
                ylim(min(ylm,ylnew));
            end
            ylabel('Cost');
        end
        title(message);
        subplot(212);   hold on;
        plot(0:it,BestCosts(1:it+1),'b','LineWidth',2);
        plot(0:it,HistoryBest(1:it+1),'r','LineWidth',1);
        if HistoryBest(it+1)>0,	set(gca,'YScale','log');	end
        xlim([0,it+1]);	grid on;
        legend({'Iteration Best','History Best'});
        xlabel('Iteration');	ylabel('Cost');
    end

    function out = output(flag)
        if ShowIterInfo
            if ShowIterInfo>1,  fprintf(message);   end
            fprintf('\n');
            disp(num2str(flag,'Exitflag: %d'));
            switch flag
                case 1,	disp('Maximum number of stall iterations reached.');
                case -1,	disp('Maximum number of iterations reached.');
            end
        end
        out.exitflag = flag;
        out.FinalIteration = it;
        out.pop = particle;
        out.BestSol = GlobalBest;
        out.BestCosts = BestCosts;
        out.iterations = 0:MaxIt;
        out.HistoryBest = HistoryBest;
        if withNeighbor
            out.NeighborNumber = [0,updateNbr];
        end
    end

end