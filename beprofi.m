function Lxy = beprofi(x, y, sigma, a)
    term_x = erf((a - x) ./ (sqrt(2) * sigma)) + erf(x ./ (sqrt(2) * sigma));
    term_y = erf((a - y) ./ (sqrt(2) * sigma)) + erf(y ./ (sqrt(2) * sigma));
    
    Lxy = (term_x .* term_y) ./ (4 * a^2);
end
