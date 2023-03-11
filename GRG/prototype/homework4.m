%%%%%%%%%%%%%% Main Entrance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% By Max Yi Ren and Emrah Bayrak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instruction: Please read through the code and fill in blanks 
% (marked by ***). Note that you need to do so for every involved
% function, i.e., m files. 

%% Optional overhead
clear; % Clear the workspace
close all; % Close all windows

%% Optimization settings
% Here we specify the objective function by giving the function handle to a
% variable, for example:
f = @(x)rosenbrock(x); % replace rosenbrock with your objective function
% In the same way, we also provide the gradient of the 
% objective:
g = @(x)rosenbrockg(x); % replace accordingly

% Specify the equality constraints
h = @(x)constraint(x);
% Provide the gradient of the constraints.
delh = @(x)constraintg(x);

% Note that explicit gradient information is only optional.
% However, providing these information to the search algorithm will save
% computational cost from finite difference calculations for them.


% Turn on or off line search. You could turn on line search once other
% parts of the program are debugged.

% opt.linesearch = ***; % false or true

% Set the tolerance to be used as a termination criterion:

% opt.eps = ***;

% Select state and decision variables and specify their indices
% E.g. if x1 and x2 are decision variables, d_id = [1,2]

% d_id = ***;
% s_id = ***;

% Set the initial guess for decision variables: (column vector)

% d0 = ***;

%% Run optimization
% Run your implementation of GRG using the gradient descent method. See
% gradient.m.

solution = gradient(f, g, h, delh, d0, d_id,s_id, opt);

%% Report
% Make sure that your solutions are saved to solution.x in your gradient.m 
% file.
report(solution,f);