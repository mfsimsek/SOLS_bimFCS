clear;
mts = 5000; %sets resolution
w = zeros(1, 30); % up to 30x30 super pixel (15x15 for bin2 acquisition)
sig = 0.38795; %PSF waist in um. Measured using sub-diffraction limit sized beads.
pix = 0.195; %camera pixel size in um (use every other w2 for bin2 acquisition).

for b = 1:30
    px = pix * b;
    stp = 0.001 * b^(1/4);
    
    coords = ((0:mts) - mts/2) * stp + px/2;
    [X, Y] = meshgrid(coords, coords);
    
    % Vectorized call to beprof
    Lmat = beprofi(X, Y, sig, px); %creates the square pixel beam profile using beprofi.m file
    
    Lmax = max(Lmat, [], 'all');
    Lmatn = Lmat / Lmax;
    
    % Use contourc to get contour matrix without plotting
    C = contourc(coords, coords, Lmatn, [1/exp(1)^2 4/exp(1)^2]);
    
    % Parse contour matrix C for the first contour level
    % C format: [level x1 x2 ...; numPoints y1 y2 ...]
    idx = 1;
    while idx < size(C, 2)
        level = C(1, idx);
        numPoints = C(2, idx);
        if abs(level - 1/exp(1)^2) < 1e-6
            x_points = C(1, idx+1 : idx+numPoints);
            y_points = C(2, idx+1 : idx+numPoints);
            break;
        end
        idx = idx + numPoints + 1;
    end
    
    xc = mean(x_points);
    yc = mean(y_points);
    dmat = sqrt((x_points - xc).^2 + (y_points - yc).^2);
    w(b) = mean(dmat);
end

w2=w.^2 %outputs w2 values calculated for the PSF and super-pixel sizes