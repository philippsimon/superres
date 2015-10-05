% Modified the code from here:
% http://stackoverflow.com/questions/2153768/draw-ellipse-and-ellipsoid-in-matlab#answer-2155162
function plotSvd(S,V)

    % range to plot over
    %------------------------------------
    N = 50;
    theta = 0:1/N:2*pi+1/N;

    % Parametric equation of the ellipse
    %----------------------------------------
    X(1,:) = S(1,1) * cos(theta); 
    X(2,:) = S(2,2) * sin(theta);

    % Coordinate transform (since your ellipse is axis aligned)
    %----------------------------------------
    X = V * X;

    % Plot
    %----------------------------------------
    plot(X(1,:),X(2,:));
end