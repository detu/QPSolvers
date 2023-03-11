% Armijo line search
function a = lineSearch(f,g,delh,x,d_id, s_id)
    
    % t = ***; % scale factor on current gradient: [0.01, 0.3]
    % b = ***; % scale factor on backtracking: [0.1, 0.8]
    % a = ***; % maximum step length
    
    % delzdeld = ***; % reduced gradient
    
    % G = ***;   % gradient vector of all variables 
    
    D = zeros(length(x),1);   % A zero direction vector for state and decision variables
    
    % Dd = ***;                  % Direction for decision variables
    % Ds = ***;                  % Direction for state variables
    D([d_id, s_id]) = [Dd;Ds];   % Save the directions to direction vector      
    
    % terminate if line search takes too long
    count = 0;
    while count<10
        % stop if condition satisfied
        
        % stop = ***;
        if stop;
            break;
        else
            % backtracking
            % a = ***;
            count = count + 1;
        end
    end