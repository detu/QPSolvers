%%%%%%%%%%%%%% Generate Report  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% By Max Yi Ren and Emrah Bayrak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function report(solution,f)
    figure; % Open an empty figure window
    grid on;
    hold on; % Hold on to the current figure
    
    % Plot the search path
    x = solution.x;
    iter = size(x,2);
    MaxMarker = 40;
    F = zeros(iter,1);
    F(1) = feval(f,x(:,1));
    Marker = 1.5^F(1);
    Marker(Marker>MaxMarker) = MaxMarker;
    plot3(x(1,1),x(2,1),x(3,1),'.b','markerSize',Marker);
    for i = 2:iter
        F(i) = feval(f,x(:,i));
        Marker = 1.5^F(i);
        % Draw lines. Type "help line" to see more drawing options.
        line([x(1,i-1),x(1,i)],[x(2,i-1),x(2,i)],[x(3,i-1),x(3,i)],'Color','b');
        plot3(x(1,i),x(2,i),x(3,i),'.b','markerSize',Marker);
    end
    plot3(x(1,i),x(2,i),x(3,i),'*r','markerSize',20);
    title(['f* = ', num2str(F(i)), ' & x* = [', num2str(x(1,i)), ', ', num2str(x(2,i)), ', ', num2str(x(3,i)), ']'])
    view(-60,20);

    % Plot the convergence
    figure;
    plot(1:iter, log(F-F(end)),'k','lineWidth',3);