function H = rosenbrockH(x)
    H = [2+800*x(1)^2+-400*(x(2)-x(1)^2),-400*x(1);
        -400*x(1),200];