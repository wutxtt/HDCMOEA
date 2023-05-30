classdef DDMOP3 < PROBLEM

    properties
        evaluated
    end

   methods
        %% Default settings of the problem
        function Setting(obj)          

            obj.M = 3;
            obj.D = obj.M+3;

            Prate = 65000;
            El = 400;
            w0 = 2 * pi * 50;
            ws = 2 * pi * 16000;
            Vdc = 800;
            Rb = El^2 / Prate;
            Iref = 141;
            K = 2;
            N = 105;
            M = K + 1;
            D = M +3;
            obj.lower = [2 * pi * Vdc / (3.2 * ws * Iref), 0.001 * Rb / w0, 0.001 / (Rb * w0), 0.001 / (Rb * w0), 0.001 / (Rb * w0) .* ones(1, M - 1)];
            obj.upper = [2 * pi * Vdc / (1.6 * ws * Iref), 0.1 * El^2 / (Prate * w0), 0.1 * El^2 / (Prate * w0), 0.05 * Prate / (El^2 * w0), 0.05 * Prate / (El^2 * w0) .* ones(1, M - 1)];
            obj.evaluated = 0;
            obj.encoding = ones(1,obj.D);
        end
%     Prate = 65000;
%     El = 400;
%     w0 = 2 * pi * 50;
%     ws = 2 * pi * 16000;
%     Vdc = 800;
%     Rb = El^2 / Prate;
%     Iref = 141;
%     K = 2;
%     N = 105;
%     M = K + 1;
%     D = M + 3;
%     lower = [2 * pi * Vdc / (3.2 * ws * Iref), 0.001 * Rb / w0, 0.001 / (Rb * w0), 0.001 / (Rb * w0), 0.001 / (Rb * w0) .* ones(1, M - 1)];
%     upper = [2 * pi * Vdc / (1.6 * ws * Iref), 0.1 * El^2 / (Prate * w0), 0.1 * El^2 / (Prate * w0), 0.05 * Prate / (El^2 * w0), 0.05 * Prate / (El^2 * w0) .* ones(1, M - 1)];
%     inverse = [2, 4, 3, 5, 1, 6];
        function PopObj = CalObj(obj,PopDec)
            if obj.evaluated <=0 
                load DDMOP3.mat initial
                inverse = [2, 4, 3, 5, 1, 6];
                PopDec(:, inverse) = initial;
                PopDec = max(min(PopDec, repmat(obj.upper, 11 * obj.D - 1, 1)), repmat(obj.lower, 11 * obj.D - 1, 1));
            end
          
            assert(obj.evaluated <= 600, 'Exceed the maximum number of function evaluations.');
            
%             PopDec = input;
            PopDec = max(min(PopDec, repmat(obj.upper, size(PopDec, 1), 1)), repmat(obj.lower, size(PopDec, 1), 1));
            PopObj = zeros(size(PopDec, 1), obj.M);
            L = PopDec(:, 1:3);
            Cf = PopDec(:, 4);
            C = PopDec(:, 5:end);
            ws = 2 * pi * 16000;
            K = 2;
            Lf = 1 ./ (repmat((1:K).^2 .* ws^2, size(C, 1), 1) .* C);

            for i = 1:obj.M - 1
                PopObj(:, i) = 20 .* log(abs(Transfer(i * ws * 1i, PopDec, obj.M - 1)));
            end

            PopObj(:, obj.M) = sum(L, 2) + sum(Lf, 2);
            pause(size(PopObj, 1) * 14 + floor(14 * 0.05 * rand))
            obj.evaluated = obj.evaluated + size(PopDec, 1);
        end
%     switch Operation
%         case 'boundary'
%             varargout = {lower, upper};
%         case 'init'
%             load DDMOP3.mat initial
%             PopDec(:, inverse) = initial;
%             PopDec = max(min(PopDec, repmat(upper, 11 * D - 1, 1)), repmat(lower, 11 * D - 1, 1));
%             evaluated = 0;
%             varargout = {PopDec};
%         case 'value'
%             assert(evaluated <= 400, 'Exceed the maximum number of function evaluations.');
% 
%             PopDec = input;
%             PopDec = max(min(PopDec, repmat(upper, size(PopDec, 1), 1)), repmat(lower, size(PopDec, 1), 1));
%             PopObj = zeros(size(PopDec, 1), M);
%             L = PopDec(:, 1:3);
%             Cf = PopDec(:, 4);
%             C = PopDec(:, 5:end);
%             Lf = 1 ./ (repmat((1:K).^2 .* ws^2, size(C, 1), 1) .* C);
% 
%             for i = 1:M - 1
%                 PopObj(:, i) = 20 .* log(abs(Transfer(i * ws * 1i, PopDec, M - 1)));
%             end
% 
%             PopObj(:, M) = sum(L, 2) + sum(Lf, 2);
%             pause(size(PopObj, 1) * 14 + floor(14 * 0.05 * rand))
%             evaluated = evaluated + size(PopDec, 1);
%             varargout = {PopObj};
    end

end

function G = Transfer(s, X, n)
    L = X(:, 1:3);
    Cf = X(:, 4);
    C = X(:, 5:end);
    ws = 2 * pi * 16000;
    Lf = 1 ./ (repmat((1:n).^2 .* ws^2, size(C, 1), 1) .* C);
    RL = 0.005;
    Glc = sum(repmat((L(:, 2) .* s .* (L(:, 3) .* Cf .* s.^2 + 1) + L(:, 3) .* s), 1, n) ./ (Lf .* s + RL + 1 ./ (C .* s)), 2);
    G = 1 ./ (L(:, 1) .* s .* Glc + L(:, 2) .* s .* (L(:, 3) .* Cf .* s.^2 + 1) + L(:, 3) .* s);
end
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP4.p
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP4.p
% b'v01.00v00.00\x00\x0bp\x1cx\xb6\xff\xb1\x00\x00\x03[\x00\x00\x03\xbb\x00\x00\x06\x82'
