classdef DDMOP7 < PROBLEM
    
%     persistent TrainIn TrainOut evaluated;
    properties
        TrainIn
        TrainOut
        evaluated
        dataNo
        nHidden
    end

    methods
        %% Default settings of the problem
        function Setting(obj)          
            obj.dataNo = 1;
            obj.nHidden = 1;

            str = {'Statlog_Australian', 'Climate', 'Statlog_German', 'Breast_cancer_Wisconsin_Diagnostic', 'Connectionist_bench_Sonar'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1). file), 'DDMOP7.mat'), 'Dataset');
            Data = Dataset.(str{obj.dataNo});
            Mean = mean(Data(:, 1:end - 1), 1);
            Std = std(Data(:, 1:end - 1), [], 1);
            Data(:, 1:end - 1) = (Data(:, 1:end - 1) - repmat(Mean, size(Data, 1), 1)) ./ repmat(Std, size(Data, 1), 1);
            Data(:, end) = Data(:, end) == Data(1, end);
            obj.TrainIn = Data(:, 1:end - 1);
            obj.TrainOut = Data(:, end);

            obj.D = (size(obj.TrainIn, 2) + 1) * obj.nHidden + (obj.nHidden + 1) * size(obj.TrainOut, 2);
%             PopDec = (rand(11 * D - 1, D) - 0.5) * 2 .* randi([0, 1], 11 * D - 1, D);
            obj.evaluated = 0;

            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = (14 + 1) * obj.nHidden + (obj.nHidden + 1) * 1; end
%             D = (14 + 1) * nHidden + (nHidden + 1) * 1;
            obj.lower = zeros(1, obj.D) - 1;
            obj.upper = zeros(1, obj.D) + 1;
%             obj.lower    = zeros(1,obj.D);
%             obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
%             g      = 100*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
%             PopObj = 0.5*repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),PopDec(:,1:obj.M-1)],2)).*[ones(size(PopDec,1),1),1-PopDec(:,obj.M-1:-1:1)];


            assert(obj.evaluated <= 600, 'Exceed the maximum number of function evaluations.');
%             PopDec = input;

            PopObj = zeros(size(PopDec, 1), 2);

            for i = 1:size(PopDec, 1)
                W1 = reshape(PopDec(i, 1:(size(obj.TrainIn, 2) + 1) * obj.nHidden), size(obj.TrainIn, 2) + 1, obj.nHidden);
                W2 = reshape(PopDec(i, (size(obj.TrainIn, 2) + 1) * obj.nHidden + 1:end), obj.nHidden + 1, size(obj.TrainOut, 2));

                [W1, W2] = Train(obj.TrainIn, obj.TrainOut, W1, W2, 1);

                Z = Predict(obj.TrainIn, W1, W2);
                PopDec(i, :) = [W1(:)', W2(:)'];

                PopObj(i, 1) = mean(PopDec(i, :) ~= 0);

                PopObj(i, 2) = mean(round(Z) ~= obj.TrainOut);
            end

            pause(size(PopObj, 1) * 6.24 + floor(6.24 * 0.05 * rand))
            obj.evaluated = obj.evaluated + size(PopDec, 1);
%             varargout = {PopObj};
        end
%         %% Generate points on the Pareto front
%         function R = GetOptimum(obj,N)
%             R = UniformPoint(N,obj.M)/2;
%         end
%         %% Generate the image of Pareto front
%         function R = GetPF(obj)
%             if obj.M == 2
%                 R = obj.GetOptimum(100);
%             elseif obj.M == 3
%                 a = linspace(0,1,10)';
%                 R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
%             else
%                 R = [];
%             end
%         end
    end
end

% function varargout = DDMOP7(varargin)
% 
%     persistent TrainIn TrainOut evaluated;
% 
%     dataNo = 1;
%     nHidden = 1;
% 
%     if nargin == 1
%         Operation = varargin{1};
%     else
%         Operation = varargin{1};
%         input = varargin{2};
%     end
% 
%     switch Operation
%         case 'boundary'
%             D = (14 + 1) * nHidden + (nHidden + 1) * 1;
%             lower = zeros(1, D) - 1;
%             upper = zeros(1, D) + 1;
%             varargout = {lower, upper};
%         case 'init'
%             str = {'Statlog_Australian', 'Climate', 'Statlog_German', 'Breast_cancer_Wisconsin_Diagnostic', 'Connectionist_bench_Sonar'};
%             CallStack = dbstack('-completenames');
%             load(fullfile(fileparts(CallStack(1). file), 'DDMOP7.mat'), 'Dataset');
%             Data = Dataset.(str{dataNo});
%             Mean = mean(Data(:, 1:end - 1), 1);
%             Std = std(Data(:, 1:end - 1), [], 1);
%             Data(:, 1:end - 1) = (Data(:, 1:end - 1) - repmat(Mean, size(Data, 1), 1)) ./ repmat(Std, size(Data, 1), 1);
%             Data(:, end) = Data(:, end) == Data(1, end);
%             TrainIn = Data(:, 1:end - 1);
%             TrainOut = Data(:, end);
% 
%             D = (size(TrainIn, 2) + 1) * nHidden + (nHidden + 1) * size(TrainOut, 2);
%             PopDec = (rand(11 * D - 1, D) - 0.5) * 2 .* randi([0, 1], 11 * D - 1, D);
%             evaluated = 0;
%             varargout = {PopDec};
%         case 'value'
%             assert(evaluated <= 600, 'Exceed the maximum number of function evaluations.');
%             PopDec = input;
% 
%             PopObj = zeros(size(PopDec, 1), 2);
% 
%             for i = 1:size(PopDec, 1)
%                 W1 = reshape(PopDec(i, 1:(size(TrainIn, 2) + 1) * nHidden), size(TrainIn, 2) + 1, nHidden);
%                 W2 = reshape(PopDec(i, (size(TrainIn, 2) + 1) * nHidden + 1:end), nHidden + 1, size(TrainOut, 2));
% 
%                 [W1, W2] = Train(TrainIn, TrainOut, W1, W2, 1);
% 
%                 Z = Predict(TrainIn, W1, W2);
%                 PopDec(i, :) = [W1(:)', W2(:)'];
% 
%                 PopObj(i, 1) = mean(PopDec(i, :) ~= 0);
% 
%                 PopObj(i, 2) = mean(round(Z) ~= TrainOut);
%             end
% 
%             pause(size(PopObj, 1) * 6.24 + floor(6.24 * 0.05 * rand))
%             evaluated = evaluated + size(PopDec, 1);
%             varargout = {PopObj};
%     end
% 
% end

function [W1, W2] = Train(X, T, W1, W2, nEpoch)

    for epoch = 1:nEpoch
        [Z, Y] = Predict(X, W1, W2);
        P = (Z - T) .* Z .* (1 - Z);
        Q = P * W2(2:end, :)' .* (1 - Y.^2);
        D1 = 0;
        D2 = 0;

        for i = 1:size(X, 1)
            D2 = D2 + [1, Y(i, :)]' * P(i, :);
            D1 = D1 + [1, X(i, :)]' * Q(i, :);
        end

        W1 = W1 - D1 / size(X, 1);
        W2 = W2 - D2 / size(X, 1);
    end

end

function [Z, Y] = Predict(X, W1, W2)
    Y = 1 - 2 ./ (1 + exp(2 * [ones(size(X, 1), 1), X] * W1));
    Z = 1 ./ (1 + exp(- [ones(size(Y, 1), 1), Y] * W2));
end
