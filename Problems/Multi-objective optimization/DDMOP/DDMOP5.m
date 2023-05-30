classdef DDMOP5 < PROBLEM

      properties
            evaluated
    end
%     lower = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0, 0, 0, 0];
%     upper = [1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.5, 0.5, 0.5, 0.5];
%     inverse = [6, 10, 9, 7, 2, 8, 5, 3, 4, 1, 11];
%     D = 11;
    
  methods
    function Setting(obj)       
        obj.lower = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0, 0, 0, 0];
        obj.upper = [1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 0.5, 0.5, 0.5, 0.5];
        inverse =  [6, 10, 9, 7, 2, 8, 5, 3, 4, 1, 11];
        obj.D = 11;
        obj.M = 3;
        load DDMOP5.mat initial
        PopDec(:, inverse) = initial;
        PopDec = PopDec .* repmat(obj.upper - obj.lower, 11 * obj.D - 1, 1) + repmat(obj.lower, 11 * obj.D - 1, 1);
        PopDec = max(min(PopDec, repmat(obj.upper, 11 * obj.D - 1, 1)), repmat(obj.lower, 11 * obj.D - 1, 1));
        obj.evaluated = 0;
        obj.encoding = ones(1,obj.D);
        
    end

    function PopObj = CalObj(obj,PopDec)
        assert(obj.evaluated <= 800, 'Exceed the maximum number of function evaluations.');
        PopDec = max(min(PopDec, repmat(obj.upper, size(PopDec, 1), 1)), repmat(obj.lower, size(PopDec, 1), 1));
        PopObj = zeros(size(PopDec, 1), 3);
        PopDec(:, 5:7) = floor((PopDec(:, 5:7) - 0.9) / 0.0125) * 0.0125 + 0.9;
        PopDec(:, 8:11) = floor((PopDec(:, 8:11)) / 0.1) * 0.1;

        for i = 1:size(PopDec, 1)
                [PopObj(i, 1), PopObj(i, 2), PopObj(i, 3)] = NR(PopDec(i, :));
        end

        obj.evaluated = obj.evaluated + size(PopDec, 1);
    end
%     switch Operation
%         case 'boundary'
%             varargout = {lower, upper};
%         case 'init'
%             load DDMOP5.mat initial
%             PopDec(:, inverse) = initial;
%             PopDec = PopDec .* repmat(upper - lower, 11 * D - 1, 1) + repmat(lower, 11 * D - 1, 1);
%             PopDec = max(min(PopDec, repmat(upper, 11 * D - 1, 1)), repmat(lower, 11 * D - 1, 1));
%             evaluated = 0;
%             varargout = {PopDec};
%         case 'value'
%             assert(evaluated <= 800, 'Exceed the maximum number of function evaluations.');
%             PopDec = input;
%             PopDec = max(min(PopDec, repmat(upper, size(input, 1), 1)), repmat(lower, size(input, 1), 1));
%             PopObj = zeros(size(PopDec, 1), 3);
%             PopDec(:, 5:7) = floor((PopDec(:, 5:7) - 0.9) / 0.0125) * 0.0125 + 0.9;
%             PopDec(:, 8:11) = floor((PopDec(:, 8:11)) / 0.1) * 0.1;
% 
%             for i = 1:size(PopDec, 1)
%                 [PopObj(i, 1), PopObj(i, 2), PopObj(i, 3)] = NR(PopDec(i, :));
%             end
% 
%             evaluated = evaluated + size(PopDec, 1);
%             varargout = {PopObj};
    end
end

function [PL, DeltaV, Vsm] = NR(PopDec)

    Data1 = [1, 3, 1.060, 0, 0, 0, 232.3859, - 16.8888, 1.06, 0;
        2, 2, 1.045, - 4.980, 21.700, 12.7, 40, 42.39645, 1.045, 0;
        3, 2, 1.010, - 12.710, 94.200, 19, 0, 23.39357, 1.01, 0;
        4, 1, 1.018, - 10.320, 47.800, - 3.9, 0, 0, 1, 0;
        5, 1, 1.020, - 8.7820, 7.600, 1.6, 0, 0, 1, 0;
        6, 2, 1.070, - 14.220, 11.200, 7.5, 0, 12.24036, 1.07, 0;
        7, 1, 1.062, - 13.360, 0, 0, 0, 0, 1, 0;
        8, 2, 1.090, - 13.360, 0, 0, 0, 17.3566, 1.09, 0;
        9, 1, 1.056, - 14.940, 29.500, 16.6, 0, 0, 1, 0.19;
        10, 1, 1.051, - 15.100, 9, 5.8, 0, 0, 1, 0;
        11, 1, 1.057, - 14.790, 3.500, 1.8, 0, 0, 1, 0;
        12, 1, 1.055, - 15.070, 6.100, 1.6, 0, 0, 1, 0;
        13, 1, 1.050, - 15.150, 13.500, 5.8, 0, 0, 1, 0;
        14, 1, 1.035, - 16.030, 14.900, 5, 0, 0, 1, 0];

    Data2 = [5, 6, 2, 0, 0.25202, 0, 0.932;
        4, 7, 2, 0, 0.20912, 0, 0.978;
        4, 9, 2, 0, 0.55618, 0, 0.969;
        1, 2, 1, 0.01938, 0.05917, 0.0528, 0;
        2, 3, 1, 0.04699, 0.19797, 0.0438, 0;
        2, 4, 1, 0.05811, 0.17632, 0.0374, 0;
        1, 5, 1, 0.05403, 0.22304, 0.0492, 0;
        2, 5, 1, 0.05695, 0.17388, 0.034, 0;
        3, 4, 1, 0.06701, 0.17103, 0.0346, 0;
        4, 5, 1, 0.01335, 0.04211, 0.0128, 0;
        7, 8, 1, 0, 0.17615, 0, 0;
        7, 9, 1, 0, 0.11001, 0, 0;
        9, 10, 1, 0.03181, 0.0845, 0, 0;
        6, 11, 1, 0.09498, 0.1989, 0, 0;
        6, 12, 1, 0.12291, 0.25581, 0, 0;
        6, 13, 1, 0.06615, 0.13027, 0, 0;
        9, 14, 1, 0.12711, 0.27038, 0, 0;
        10, 11, 1, 0.08205, 0.19207, 0, 0;
        12, 13, 1, 0.22092, 0.19988, 0, 0;
        13, 14, 1, 0.17093, 0.34802, 0, 0];

    Data1([2, 3, 6, 8], 9) = PopDec(1:4);
    Data2([1, 2, 3], 7) = PopDec(5:7);
    Data1([9, 10, 13, 14], 10) = PopDec(8:11);
    baseMVA = 100;
    Bus = Data1(:, 1);
    Vtype = Data1(:, 2);
    Pload = Data1(:, 5);
    Qload = Data1(:, 6);
    Pgen = Data1(:, 7);
    Qgen = Data1(:, 8);
    Vset = Data1(:, 9);
    Qsh = Data1(:, 10);
    II = Data2(:, 1);
    JJ = Data2(:, 2);
    Ltype = Data2(:, 3);
    R = Data2(:, 4);
    X = Data2(:, 5);
    B = Data2(:, 6) / 2;
    K = Data2(:, 7);

    y1 = zeros(14);
    y2 = zeros(14);
    y3 = zeros(14);
    lin = length(II);

    for x = 1:lin

        switch Ltype(x)
            case 1
                y1(II(x), JJ(x)) = 1 / (R(x) + 1i * X(x));
                y1(JJ(x), II(x)) = y1(II(x), JJ(x));
                y3(II(x), JJ(x)) = 1i * B(x);
                y3(JJ(x), II(x)) = 1i * B(x);
            case 2
                y1(II(x), JJ(x)) = 1 / ((R(x) + 1i * X(x)) * K(x));
                y1(JJ(x), II(x)) = y1(II(x), JJ(x));
                y2(II(x), JJ(x)) = (1 - K(x)) / ((R(x) + 1i * X(x)) * K(x)^2);
                y2(JJ(x), II(x)) = (K(x) - 1) / ((R(x) + 1i * X(x)) * K(x));
        end

    end

    Y = zeros(14);

    for x = 1:14
        Y(x, x) = sum(y1(x, :)) + sum(y2(x, :)) + sum(y3(x, :)) + 1i * Qsh(x);
    end

    Y = Y - y1;
    G = real(Y);
    B = imag(Y);

    U = Vset;
    e = real(U);
    f = imag(U);

    Ps = zeros(1, 14);
    Qs = zeros(1, 14);
    D = ones(26, 1);

    for x = 1:14
        Ps(x) = (Pgen(x) - Pload(x)) / baseMVA;
        Qs(x) = (Qgen(x) - Qload(x)) / baseMVA;
    end

    N = 0;
    Jacbi = zeros(26);

    while max(abs(D)) > 0.000001

        for x = 2:14

            switch Vtype(x)
                case 1
                    D(2 * x - 3) = Ps(x) - e(x) * (G(x, :) * e - B(x, :) * f) - f(x) * (G(x, :) * f + B(x, :) * e);
                    D(2 * x - 2) = Qs(x) - f(x) * (G(x, :) * e - B(x, :) * f) + e(x) * (G(x, :) * f + B(x, :) * e);
                case 2
                    D(2 * x - 3) = Ps(x) - e(x) * (G(x, :) * e - B(x, :) * f) - f(x) * (G(x, :) * f + B(x, :) * e);
                    D(2 * x - 2) = Vset(x) .* Vset(x) - (e(x).^2 + f(x).^2);
            end

        end

        for I = 2:14

            for J = 2:14

                if I ~= J
                    Jacbi((2 * I - 3), (2 * J - 3)) = B(I, J) * e(I) - G(I, J) * f(I);
                    Jacbi((2 * I - 3), (2 * J - 2)) =- (G(I, J) * e(I) + B(I, J) * f(I));

                    switch Vtype(I)
                        case 1
                            Jacbi((2 * I - 2), (2 * J - 3)) = G(I, J) * e(I) + B(I, J) * f(I);
                            Jacbi((2 * I - 2), (2 * J - 2)) = B(I, J) * e(I) - G(I, J) * f(I);
                        case 2
                            Jacbi((2 * I - 2), (2 * J - 3)) = 0;
                            Jacbi((2 * I - 2), (2 * J - 2)) = 0;
                    end

                else
                    Jacbi(2 * I - 3, 2 * J - 3) =- (G(I, :) * f + B(I, :) * e) + B(I, I) * e(I) - G(I, I) * f(I);
                    Jacbi(2 * I - 3, 2 * J - 2) =- (G(I, :) * e - B(I, :) * f) - G(I, I) * e(I) - B(I, I) * f(I);

                    switch Vtype(I)
                        case 1
                            Jacbi(2 * I - 2, 2 * J - 3) =- (G(I, :) * e - B(I, :) * f) + G(I, I) * e(I) + B(I, I) * f(I);
                            Jacbi(2 * I - 2, 2 * J - 2) = (G(I, :) * f + B(I, :) * e) + B(I, I) * e(I) - G(I, I) * f(I);
                        case 2
                            Jacbi(2 * I - 2, 2 * J - 3) =- 2 * f(I);
                            Jacbi(2 * I - 2, 2 * J - 2) =- 2 * e(I);
                    end

                end

            end

        end

        Deta =- inv(Jacbi) * D;

        for x = 2:14
            f(x) = f(x) + Deta((2 * x - 3), 1);
            e(x) = e(x) + Deta((2 * x - 2), 1);
        end

        U = e + 1i * f;
        N = N + 1;
    end

    N = N - 1;
    S0 = U(1) * (conj(Y(1, :)) * conj(U));
    S1 = zeros(20, 1);
    S2 = zeros(20, 1);

    for x = 1:20
        S1(x) = U(II(x)) * (conj(U(II(x))) * (conj(y2(II(x))) + y3(II(x)) + 1i * Qsh(II(x))) +...
        (conj(U(II(x))) - conj(U(JJ(x)))) * conj(y1(II(x), JJ(x))));
        S2(x) = U(JJ(x)) * (conj(U(JJ(x))) * (conj(y2(JJ(x))) + y3(JJ(x)) + 1i * Qsh(JJ(x))) +...
        (conj(U(JJ(x))) - conj(U(II(x)))) * conj(y1(II(x), JJ(x))));
    end

    detaS = S1 + S2;
    Vabs = abs(U);
    Angle = atan(f ./ e) * 180 / pi;
    PL = sum(real(detaS));
    Vsm = 1 / min(svd(Jacbi));
    phi = (abs(Vabs(Data1(:, 2) == 1) - 1) - 0.05) ./ Vabs(Data1(:, 2) == 1);
    DeltaV = sum(phi(phi > 0));

end
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP6.p
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP6.p
% b'v01.00v00.00\x00\x01\x10\x1c\xc0\xaf\xbf\xb1\x00\x00\x01\xa2\x00\x00\x02\x9b\x00\x00\x03\xf5'
