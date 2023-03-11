function s = solveh(x, h, delh, s_id)

eps = 1e-8;  % Set a tolerance for convergence 
% s = ***;    % Save the current state variables

iter = 0;       % Set initial iteration to 0
% Set termination criterion
% hnorm = ***; % norm of the constraint vector

    while(hnorm > eps)
        iter = iter+1; % Increase iteration by 1
        
        % dhds = ***; % current dh/ds
        
        % Modify dh/ds when it is singula
        %%% KEEP THIS %%%
        dhds_inv = correctH(dhds);
        %%%%%%%%%%%%%%%%%
        
        % s = ***; % Use modified dh/ds to calculate new s
        x(s_id) = s;        % Save new s to the current solution
        
        % hnorm = ***; % Update termination critetion
    end
end