
%{
Variables and what they mean
i: The current lag/ x-axis position
- Represents the actual lag time in frames.

n: the "power" of the bin/ the multiplier
- controls the exponent for logarithmic binning
- as n increases, the bin gets wider

k: the result index
- a counter for the output array

j: counter for averaging
- exists within the binning loop for averaging values within bins.

%}

function [curve_log,lagtimeAC] = lintolog(curve)
    curve_log = nan(150, 1); % Preallocate curve_log
    lagtimeAC = nan(150, 1); % Preallocate lagtimeAC

    k = 1; % Initialize index for lagtimeAC
    for i = 1:10
        curve_log(i) = curve(round(i + 0.5 * (length(curve) - 1))); % Round to nearest integer
        lagtimeAC(k) = i;
        k = k + 1;
    end

    n = 1; % Initialize n

    while i <= 0.5 * (length(curve) - 1)-2^n
        t = 10; % Initialize t
        while t > 0
            % Calculate the start and end indices for the current bin
            start_idx = floor(i + 0.5 * (length(curve) - 1));
            end_idx = start_idx + (2^n - 1);

            % Check if the ENTIRE bin fits within the array
            if end_idx <= length(curve)
                curve_log(k) = 0; % Reset curve_log for current k
                for j = 0:(2^n - 1)
                    % Use the start_idx we already calculated
                    curve_log(k) = curve_log(k) + curve(start_idx + j);
                end
                curve_log(k) = curve_log(k) / (2^n); % Average the values
                lagtimeAC(k) = i + 0.5 * (2^n);
                i = i + 2^n; % Increment i
                k = k + 1; % Increment k
                t = t - 1; % Decrement t
            else
                % If the bin would go out of bounds, stop binning immediately
                t = 0;
                break;
            end   %terminate if statement block
        end
        n = n + 1; % Increment n
    end
    lagtimeAC(1)=nan;
    curve_log(1)=nan;
    %curve_log(k:end) = []; % Remove unused preallocated entries
    %lagtimeAC(k:end) = []; % Remove unused preallocated entries
end
