function varargout = DDMOP1(varargin)
    persistent evaluated;
    M = 9;
    D = 11;
    lower = [0.5, 0.45, 0.5, 0.5, 0.875, 0.4, 0.4, - 1.655, - 1.808, - 2, - 2];
    upper = [1.5, 1.35, 1.5, 1.5, 2.625, 1.2, 1.2, 2.345, 2.192, 2, 2];
    inverse = [6, 3, 11, 7, 8, 5, 1, 2, 4, 9, 10];

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
            load 'DDMOP1.mat' intial
            PopDec(:, inverse) = intial;
            PopDec = PopDec .* repmat(upper - lower, 11 * D - 1, 1) + repmat(lower, 11 * D - 1, 1);
            PopDec = max(min(PopDec, repmat(upper, 11 * D - 1, 1)), repmat(lower, 11 * D - 1, 1));
            varargout = {PopDec};
            evaluated = 0;
        case 'value'

            assert(evaluated <= 400, 'Exceed the maximum number of function evaluations.');

            PopDec = input;
            PopDec = max(min(PopDec, repmat(upper, size(input, 1), 1)), repmat(lower, size(input, 1), 1));
            PopObj = zeros(size(PopDec, 1), M);
            PopObj(:, 1) = 1.98 + 4.9 .* PopDec(:, 1) + 6.67 .* PopDec(:, 2) + 6.98 .* PopDec(:, 3)
            + 4.01 .* PopDec(:, 4) + 1.78 .* PopDec(:, 5) + 0.00001 .* PopDec(:, 5) + 2.73 .* PopDec(:, 2);
            PopObj(:, 2) = 1.16 - 0.3717 .* PopDec(:, 2) .* PopDec(:, 4) - 0.00931 .* PopDec(:, 2) .* PopDec(:, 10)
            - 0.484 .* PopDec(:, 3) .* PopDec(:, 9) + 0.01343 .* PopDec(:, 6) .* PopDec(:, 10);
            PopObj(:, 3) = 0.261 - 0.0159 .* PopDec(:, 1) .* PopDec(:, 2) - 0.188 .* PopDec(:, 1) .* PopDec(:, 8)
            - 0.019 .* PopDec(:, 2) .* PopDec(:, 7) + 0.0144 .* PopDec(:, 3) .* PopDec(:, 5) + 0.87570 .* PopDec(:, 5) .* PopDec(:, 10)
            + 0.08045 .* PopDec(:, 6) .* PopDec(:, 9) + 0.00139 .* PopDec(:, 8) .* PopDec(:, 11) + 0.00001575 .* PopDec(:, 10) .* PopDec(:, 11);
            PopObj(:, 4) = 0.214 + 0.00817 .* PopDec(:, 5) - 0.131 .* PopDec(:, 1) .* PopDec(:, 8) - 0.0704 .* PopDec(:, 1) .* PopDec(:, 9)
            + 0.03099 .* PopDec(:, 2) .* PopDec(:, 6) - 0.018 .* PopDec(:, 2) .* PopDec(:, 7) + 0.0208 .* PopDec(:, 3) .* PopDec(:, 8)
            + 0.121 .* PopDec(:, 3) .* PopDec(:, 9) - 0.00364 .* PopDec(:, 5) .* PopDec(:, 6) + 0.0007715 .* PopDec(:, 5) .* PopDec(:, 10)
            - 0.0005354 .* PopDec(:, 6) .* PopDec(:, 10) + 0.00121 .* PopDec(:, 8) .* PopDec(:, 11) + 0.00184 .* PopDec(:, 9) .* PopDec(:, 10)
            - 0.018 .* PopDec(:, 2).^2;
            PopObj(:, 5) = 0.74 - 0.61 .* PopDec(:, 2) - 0.163 .* PopDec(:, 3) .* PopDec(:, 8)
            + 0.001232 .* PopDec(:, 3) .* PopDec(:, 10) - 0.166 .* PopDec(:, 7) .* PopDec(:, 9)
            + 0.227 .* PopDec(:, 2).^2;
            PopObj(:, 6) = 109.2 - 9.9 .* PopDec(:, 2) + 6.768 .* PopDec(:, 3) + 0.1792 .* PopDec(:, 10)
            - 9.256 .* PopDec(:, 1) .* PopDec(:, 2) - 12.9 .* PopDec(:, 1) .* PopDec(:, 8) - 11 .* PopDec(:, 2) .* PopDec(:, 8)
            + 0.1107 .* PopDec(:, 3) .* PopDec(:, 10) + 0.0207 .* PopDec(:, 5) .* PopDec(:, 10) + 6.63 .* PopDec(:, 6) .* PopDec(:, 9)
            - 17.75 .* PopDec(:, 7) .* PopDec(:, 8) + 22 .* PopDec(:, 8) .* PopDec(:, 9) + 0.32 .* PopDec(:, 9) .* PopDec(:, 10);
            PopObj(:, 7) = 4.72 - 0.5 .* PopDec(:, 4) - 0.19 .* PopDec(:, 2) .* PopDec(:, 3) - 0.0122 .* PopDec(:, 4) .* PopDec(:, 10)
            + 0.009325 .* PopDec(:, 6) .* PopDec(:, 10) + 0.000191 .* PopDec(:, 11).^2;
            PopObj(:, 8) = 10.58 - 0.674 .* PopDec(:, 1) .* PopDec(:, 2) - 1.95 .* PopDec(:, 2) .* PopDec(:, 8) + 0.02054 .* PopDec(:, 3) .* PopDec(:, 10)
            - 0.0198 .* PopDec(:, 4) .* PopDec(:, 10) + 0.028 .* PopDec(:, 6) .* PopDec(:, 10);
            PopObj(:, 9) = 16.45 - 0.489 .* PopDec(:, 3) .* PopDec(:, 7) - 0.843 .* PopDec(:, 5) .* PopDec(:, 6) + 0.0432 .* PopDec(:, 9) .* PopDec(:, 10);

            evaluated = evaluated + size(PopDec, 1);
            pause(size(PopObj, 1) * 16 + floor(16 * 0.05 * rand))
            varargout = {PopObj};
    end

end
E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP2.p
E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP2.p
b'v01.00v00.00\x00\x0e\xb0\x1cC\xc6_\xb1\x00\x00\x02\x8d\x00\x00\x03\t\x00\x00\x05\x84'
