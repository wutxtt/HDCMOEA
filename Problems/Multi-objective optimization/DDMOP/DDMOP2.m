classdef DDMOP2 < PROBLEM
  
    properties
        evaluated
    end
    methods
        %% Default settings of the problem
        function Setting(obj)         
            obj.M = 3;
            obj.D = 5;
            obj.lower = ones(1,obj.D);
            obj.upper = 3.*ones(1,obj.D);
            obj.evaluated = 0;
            obj.encoding = ones(1,obj.D);
        end
%     persistent evaluated;
%     M = 3;
%     D = 5;
%     lower = ones(1, D);
%     upper = 3 .* ones(1, D);
%     inverse = [3, 5, 1, 2, 4];
% 
%     if nargin == 1
%         Operation = varargin{1};
%     else
%         Operation = varargin{1};
%         input = varargin{2};
%     end
        function PopObj = CalObj(obj,PopDec)
%             if obj.evaluated <=0 
%                 load DDMOP2.mat initial
%                 inverse = [3, 5, 1, 2, 4];
%                 PopDec(:, inverse) = initial;
%                 PopDec = PopDec .* repmat(obj.upper - obj.lower, 11 * obj.D - 1, 1) + repmat(obj.lower, 11 * obj.D - 1, 1);
%                 PopDec = max(min(PopDec, repmat(obj.upper, 11 * obj.D - 1, 1)), repmat(obj.lower, 11 * obj.D - 1, 1));
%                 
%             end
%             assert(obj.evaluated <= 1000, 'Exceed the maximum number of function evaluations.');
            PopDec = max(min(PopDec, repmat(obj.upper, size(PopDec, 1), 1)), repmat(obj.lower, size(PopDec, 1), 1));
            PopObj = zeros(size(PopDec, 1), obj.M);
            PopObj(:, 1) = 1640.2823 + 2.3573285 .* PopDec(:, 1) + 2.3220035 .* PopDec(:, 2)...
            + 4.5688768 .* PopDec(:, 3) + 7.7213633 .* PopDec(:, 4) + 4.4559504 .* PopDec(:, 5);
            PopObj(:, 2) = 6.5856 + 1.15 .* PopDec(:, 1) - 1.0427 .* PopDec(:, 2)...
            + 0.9738 .* PopDec(:, 3) + 0.8364 .* PopDec(:, 4) - 0.3695 .* PopDec(:, 1) .* PopDec(:, 4)...
            + 0.0861 .* PopDec(:, 1) .* PopDec(:, 5) + 0.3628 .* PopDec(:, 2) .* PopDec(:, 4) - 0.1106 .* PopDec(:, 1).^2 ...
            - 0.3437 .* PopDec(:, 3).^2 + 0.1764 .* PopDec(:, 4).^2;
            PopObj(:, 3) =- 0.0551 + 0.0181 .* PopDec(:, 1) + 0.1024 .* PopDec(:, 2)...
            + 0.0421 .* PopDec(:, 3) - 0.0073 .* PopDec(:, 1) .* PopDec(:, 2) + 0.024 .* PopDec(:, 2) .* PopDec(:, 3)...
            - 0.0118 .* PopDec(:, 2) .* PopDec(:, 4) - 0.0204 .* PopDec(:, 3) .* PopDec(:, 4) - 0.008 .* PopDec(:, 3) .* PopDec(:, 5)...
            - 0.0241 .* PopDec(:, 2).^2 + 0.0109 .* PopDec(:, 4).^2;
            pause(size(PopObj, 1) * 19 + floor(19 * 0.05 * rand))

            obj.evaluated = obj.evaluated + size(PopDec, 1);
        end
        function R = GetOptimum(obj,N)
            R = [1689.74 5.59567 0.055533;
                    1696.3 5.55174 0.066563;
                    1675.58 4.47847 0.262609;
                    1690.62 6.357 0.050618;
                    1669.14 6.67743 0.062413;
                    1697.93 5.64378 0.055403;
                    1695 5.72354 0.039277;
                    1679.62 4.88909 0.212811;
                    1665.69 6.42572 0.068167;
                    1672.03 6.51818 0.060617;
                    1688.4 5.07193 0.165519;
                    1687.36 6.4424 0.055686;
                    1678.59 5.67282 0.061767;
                    1681.79 4.63313 0.227072;
                    1666.89 5.30448 0.100443;
                    1661.71 6.63164 0.070684;
                    1670.62 6.80367 0.05471;
                    1693.31 5.68448 0.042852;
                    1667.9 6.75112 0.059584;
                    1685.35 6.1099 0.059836;
                    1682.55 4.6518 0.221514;
                    1673.2 6.46793 0.061608;
                    1691.55 5.65053 0.052289;
                    1689.11 6.4764 0.052537;
                    1664.01 6.58728 0.070601;
                    1676.82 5.72882 0.074664;
                    1675.28 4.55464 0.256598;
                    1694.99 5.72333 0.039335
                     ];
        end
%     switch Operation
%         case 'boundary'
%             varargout = {lower, upper};
%         case 'init'
%             load DDMOP2.mat initial
%             PopDec(:, inverse) = initial;
%             PopDec = PopDec .* repmat(upper - lower, 11 * D - 1, 1) + repmat(lower, 11 * D - 1, 1);
%             PopDec = max(min(PopDec, repmat(upper, 11 * D - 1, 1)), repmat(lower, 11 * D - 1, 1));
%             evaluated = 0;
%             varargout = {PopDec};
%         case 'value'
%             assert(evaluated <= 300, 'Exceed the maximum number of function evaluations.');
% 
%             PopDec = input;
%             PopDec = max(min(PopDec, repmat(upper, size(input, 1), 1)), repmat(lower, size(input, 1), 1));
%             PopObj = zeros(size(PopDec, 1), M);
%             PopObj(:, 1) = 1640.2823 + 2.3573285 .* PopDec(:, 1) + 2.3220035 .* PopDec(:, 2)
%             + 4.5688768 .* PopDec(:, 3) + 7.7213633 .* PopDec(:, 4) + 4.4559504 .* PopDec(:, 5);
%             PopObj(:, 2) = 6.5856 + 1.15 .* PopDec(:, 1) - 1.0427 .* PopDec(:, 2)
%             + 0.9738 .* PopDec(:, 3) + 0.8364 .* PopDec(:, 4) - 0.3695 .* PopDec(:, 1) .* PopDec(:, 4)
%             + 0.0861 .* PopDec(:, 1) .* PopDec(:, 5) + 0.3628 .* PopDec(:, 2) .* PopDec(:, 4) - 0.1106 .* PopDec(:, 1).^2
%             - 0.3437 .* PopDec(:, 3).^2 + 0.1764 .* PopDec(:, 4).^2;
%             PopObj(:, 3) =- 0.0551 + 0.0181 .* PopDec(:, 1) + 0.1024 .* PopDec(:, 2)
%             + 0.0421 .* PopDec(:, 3) - 0.0073 .* PopDec(:, 1) .* PopDec(:, 2) + 0.024 .* PopDec(:, 2) .* PopDec(:, 3)
%             - 0.0118 .* PopDec(:, 2) .* PopDec(:, 4) - 0.0204 .* PopDec(:, 3) .* PopDec(:, 4) - 0.008 .* PopDec(:, 3) .* PopDec(:, 5)
%             - 0.0241 .* PopDec(:, 2).^2 + 0.0109 .* PopDec(:, 4).^2;
%             pause(size(PopObj, 1) * 19 + floor(19 * 0.05 * rand))
% 
%             evaluated = evaluated + size(PopDec, 1);
%             varargout = {PopObj};
    
    end
end
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP3.p
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP3.p
% b'v01.00v00.00\x00\x0f0\x1cB&\x7f\xb1\x00\x00\x03.\x00\x00\x03\x8b\x00\x00\x06+'
