function cmap = redblue256()
    % Red-White-Blue colormap smooth a 256 livelli
    m = 256;
    bottom = [0 0 0.5];
    botmiddle = [0 0.5 1];
    middle = [1 1 1];
    topmiddle = [1 0 0];
    top = [0.5 0 0];

    % interpolate 5 control points to 256
    x = [0 0.25 0.5 0.75 1];
    r = interp1(x, [bottom(1) botmiddle(1) middle(1) topmiddle(1) top(1)], ...
        linspace(0,1,m));
    g = interp1(x, [bottom(2) botmiddle(2) middle(2) topmiddle(2) top(2)], ...
        linspace(0,1,m));
    b = interp1(x, [bottom(3) botmiddle(3) middle(3) topmiddle(3) top(3)], ...
        linspace(0,1,m));

    cmap = [r(:) g(:) b(:)];
end