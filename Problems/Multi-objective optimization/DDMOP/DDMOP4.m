function varargout = DDMOP4(varargin)
    persistent evaluated;
    Prate = 65000;
    El = 400;
    w0 = 2 * pi * 50;
    ws = 2 * pi * 16000;
    Vdc = 800;
    Rb = El^2 / Prate;
    Iref = 141;
    K = 9;
    N = 256;
    M = K + 1;
    D = M + 3;
    lower = [2 * pi * Vdc / (3.2 * ws * Iref), 0.001 * Rb / w0, 0.001 / (Rb * w0), 0.001 / (Rb * w0), 0.001 / (Rb * w0) .* ones(1, M - 1)];
    upper = [2 * pi * Vdc / (1.6 * ws * Iref), 0.1 * El^2 / (Prate * w0), 0.1 * El^2 / (Prate * w0), 0.05 * Prate / (El^2 * w0), 0.05 * Prate / (El^2 * w0) .* ones(1, M - 1)];
    inverse = [12, 5, 3, 13, 10, 8, 2, 6, 9, 4, 1, 7, 11];

    if nargin == 1
        Operation = varargin{1};
    else
        Operation = varargin{1};
        input = varargin{2};
    end

    switch Operation
        case 'boundary'
            varargout = {lower, upper};
        case 'init'
            load DDMOP4.mat initial
            PopDec(:, inverse) = initial;
            PopDec = PopDec .* repmat(upper - lower, 11 * D - 1, 1) + repmat(lower, 11 * D - 1, 1);
            PopDec = max(min(PopDec, repmat(upper, 11 * D - 1, 1)), repmat(lower, 11 * D - 1, 1));
            evaluated = 0;
            varargout = {PopDec};
        case 'value'
            assert(evaluated <= 600, 'Exceed the maximum number of function evaluations.');

            PopDec = input;
            PopDec = max(min(PopDec, repmat(upper, size(input, 1), 1)), repmat(lower, size(input, 1), 1));
            PopObj = zeros(size(PopDec, 1), M);
            L = PopDec(:, 1:3);
            Cf = PopDec(:, 4);
            C = PopDec(:, 5:eend);
            Lf = 1 ./ (repmat((1:K).^2 .* ws^2, size(C, 1), 1) .* C);

            for i = 1:M - 1
                PopObj(:, i) = 20 .* log(abs(Transfer(i * ws * 1i, PopDec, M - 1)));
            end

            PopObj(:, M) = sum(L, 2) + sum(Lf, 2);
            pause(size(PopObj, 1) * 8.5 + floor(8.5 * 0.05 * rand))
            evaluated = evaluated + size(PopDec, 1);
            varargout = {PopObj};
    end

end

function G = Transfer(s, X, n)
    L = X(:, 1:3);
    Cf = X(:, 4);
    C = X(:, 5:eend);
    ws = 2 * pi * 16000;
    Lf = 1 ./ (repmat((1:n).^2 .* ws^2, size(C, 1), 1) .* C);
    RL = 0.005;
    Glc = sum(repmat((L(:, 2) .* s .* (L(:, 3) .* Cf .* s.^2 + 1) + L(:, 3) .* s), 1, n) ./ (Lf .* s + RL + 1 ./ (C .* s)), 2);
    G = 1 ./ (L(:, 1) .* s .* Glc + L(:, 2) .* s .* (L(:, 3) .* Cf .* s.^2 + 1) + L(:, 3) .* s);
end
E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP5.p
E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP5.p
b'v01.00v00.00\x00\x04\xa0\x1c\x9eC\x1f\xb1\x00\x00\x0c!\x00\x00\t\x9d\x00\x00\x15\xc2'
