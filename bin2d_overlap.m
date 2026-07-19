
function [B, rows, cols] = bin2d_overlap(A, n)
% BIN2D_OVERLAP Average n-by-n blocks with 50% overlap.
%   [B, rows, cols] = bin2d_overlap(A, n)
%   A      - input 3-D matrix
%   n      - block size (positive integer)
%   B      - matrix of block averages
%   rows   - row indices of the top-left corner of each block
%   cols   - column indices of the top-left corner of each block

step = ceil(n/2);              % 50% overlap; for even n this is exact

H = size(A,1);
W = size(A,2);
F = size(A,3);
row_starts = 1:step:(H - n + 1);
col_starts = 1:step:(W - n + 1);

%{
 Optionally include last partial blocks to cover edges:
if isempty(row_starts) || row_starts(end) ~= (H - n + 1)
    row_starts = [row_starts, max(1, H - n + 1)];
end
if isempty(col_starts) || col_starts(end) ~= (W - n + 1)
    col_starts = [col_starts, max(1, W - n + 1)];
end
%}

R = numel(row_starts);
C = numel(col_starts);
B = nan(R, C, F);

for i = 1:R
    r = row_starts(i);
    for j = 1:C
        c = col_starts(j);
        block = A(r:(r+n-1), c:(c+n-1),:);
        %for t=1:F
        %B(i,j,t) = mean2(block(:,:,t));
        %end
        B(i,j,:) = reshape(mean(reshape(block, [], size(block,3)), 1), [1 1 size(block,3)]);
    end
end

rows = row_starts;
cols = col_starts;
end