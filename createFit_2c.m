function [p_fit, y_fit] = createFit_2c(fitx, fity, area, sig , px, minbin)
    
   a=1000*px*((1:numel(fitx))+minbin-1);
   ari=area.*10^6;

    % Initial guess: [N1, Df, Ds_1, N2_1, base_1, Ds_2, N2_2, base_2....]
    p0 = [2e2 1e7];
    % Lower and upper bounds
    lb = [0 1e2];
    ub = [2e10 2e8];
    for k=1:numel(fitx)
        p0app=[1e1, 1000/k, 0];
        lbapp=[0,3,-200];
        ubapp=[2e7, 2000, 200];
        p0=[p0, p0app];
        lb=[lb, lbapp];
        ub=[ub, ubapp];
    end

    p_binf=zeros(4,8);
    pfi=[p0(1:5),a(1),sig,ari(1)];
    lbi=[lb(1:5),a(1),sig,ari(1)];
    ubi=[ub(1:5),a(1),sig,ari(1)];
    fitxi=fitx{1};
    fityi=fity{1};
    p_binf(1,:)=lsqcurvefit(@model_fun2, pfi, fitxi, fityi, lbi, ubi);

    for i = 2:3
        pfi=[pfi(1:5),a(i),sig,ari(i)];
        lbi=[lbi(1:5),a(i),sig,ari(i)];
        ubi=[ubi(1:5),a(i),sig,ari(i)];
        fitxi=fitx{i};
        fityi=fity{i};
        p_binf(i,:)=lsqcurvefit(@model_fun2, pfi, fitxi, fityi, lbi, ubi);
        pfi(1:5)=mean(p_binf(1:i,1:5),1);           
    end

    p0(1:5)=pfi(1:5);
    lb(1)=pfi(1);
    ub(1)=pfi(1);
    lb(2)=pfi(2);
    ub(2)=pfi(2);

    % Fit using lsqcurvefit
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 2e7, 'MaxIterations', 2e7);

    % Flatten all xData into one vector for lsqcurvefit
    x_all = {fitx, ari, a, sig, fity};
    ydata=cell2mat(x_all{5});
    ydata=ydata(:);

    p_fit = lsqcurvefit(@stacked_model, p0, x_all, ydata, lb, ub, opts);

    for i = 1:numel(fitx)
        fitxi=fitx{i};
        y_fit{i} = model_fun(fitxi, p_fit(3*i), p_fit(2), p_fit(1), ...
                          p_fit(3*i+1), p_fit(3*i+2), a(i), sig, ari(i));
    end    
end

% --- Model for one dataset ---
function y = model_fun(x, Ds, Df, N1, N2, base, a, sigma, are)
    term1 = ((2*sqrt(sigma^2 + Ds.*x) .* ...
        (exp(-0.25*a^2 ./ (sigma^2 + Ds.*x)) - 1) ./ (a*sqrt(pi)) + ...
        erf(0.5*a ./ sqrt(sigma^2 + Ds.*x))).^2 ./ ...
        ( (2*sigma*(exp(-0.25*a^2/sigma^2) - 1) / (a*sqrt(pi)) + ...
        erf(0.5*a/sigma))^2 ));

    term2 = (2*sqrt(sigma^2 + Df.*x) .* ...
        (exp(-0.25*a^2 ./ (sigma^2 + Df.*x)) - 1) ./ (a*sqrt(pi)) + ...
        erf(0.5*a ./ sqrt(sigma^2 + Df.*x))).^2 ./ ...
        ( (2*sigma*(exp(-0.25*a^2/sigma^2) - 1) / (a*sqrt(pi)) + ...
        erf(0.5*a/sigma))^2 );

    y = base + (N1/are).*term1 + (N2-base).*term2;
end

% --- Stacked model for all datasets ---
function y_all = stacked_model(p, x_all)
    fitx = x_all{1};
    ari = x_all{2};
    a = x_all{3};
    sigma = x_all{4};
    %fity = x_all{4};

    Df = p(2);
    N1 = p(1);
    y_all = [];

    % Loop over datasets
    for i = 1:numel(fitx)
        Ds = p((i-1)*3 + 3);
        N2 = p((i-1)*3 + 4);
        base = p((i-1)*3 + 5);
        xdata=fitx{i};
        xdata(isnan(xdata))=[];
        y_model = model_fun(xdata, Ds, Df, N1, N2, base, a(i), sigma,ari(i));
        y_all = [y_all; y_model];
    end
end


% --- Model for one dataset ---
function y = model_fun2(p,x)
    Ds=p(3);
    Df=p(2);
    N1=p(1);
    N2=p(4);
    base=p(5);
    a=p(6);
    sigma=p(7);
    are=p(8);

    term1 = (2*sqrt(sigma^2 + Ds.*x) .* ...
        (exp(-0.25*a^2 ./ (sigma^2 + Ds.*x)) - 1) ./ (a*sqrt(pi)) + ...
        erf(0.5*a ./ sqrt(sigma^2 + Ds.*x))).^2 ./ ...
        ( (2*sigma*(exp(-0.25*a^2/sigma^2) - 1) / (a*sqrt(pi)) + ...
        erf(0.5*a/sigma))^2 );

    term2 = (2*sqrt(sigma^2 + Df.*x) .* ...
        (exp(-0.25*a^2 ./ (sigma^2 + Df.*x)) - 1) ./ (a*sqrt(pi)) + ...
        erf(0.5*a ./ sqrt(sigma^2 + Df.*x))).^2 ./ ...
        ( (2*sigma*(exp(-0.25*a^2/sigma^2) - 1) / (a*sqrt(pi)) + ...
        erf(0.5*a/sigma))^2 );

    y = base + (N1/are).*term1 + (N2-base).*term2;
end

