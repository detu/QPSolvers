function Hinv = correctH(H)
    % Fix singular H matrix
    I = eye(size(H));
    lambda = 0.01;
    Hinv = (H'*H+lambda*I)\H'; %use Levenberg-Marquardt when singular
end