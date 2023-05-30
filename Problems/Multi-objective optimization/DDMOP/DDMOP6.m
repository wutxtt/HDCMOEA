classdef DDMOP6 < PROBLEM

    properties
            Yield
            Risk
            evaluated
            dataNo
    end


   methods
        %% Default settings of the problem
        function Setting(obj)       


            obj.dataNo = 1;
            obj.D = 10;
            obj.lower = zeros(1, obj.D) - 1;
            obj.upper = zeros(1, obj.D) + 1;
            str = {'Data1', 'Data2', 'Data3', 'Data4', 'Data5'};
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1). file), 'DDMOP6.mat'), 'Dataset');
            Data = Dataset.(str{obj.dataNo});
            obj.Yield = log(Data(:, 2:end)) - log(Data(:, 1:end - 1));
            obj.Risk = cov(obj.Yield');
            obj.evaluated = 0;
            
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = size(obj.Yield, 1); end
            obj.D = size(obj.Yield, 1);
            obj.encoding = ones(1,obj.D);
        end
        
        function PopObj = CalObj(obj,PopDec)
            assert(obj.evaluated <= 300, 'Exceed the maximum number of function evaluations.');

            PopDec = PopDec ./ repmat(max(sum(abs(PopDec), 2), 1), 1, size(PopDec, 2));

            PopObj = zeros(size(PopDec, 1), 2);

            for i = 1:size(PopDec, 1)

                PopObj(i, 1) = PopDec(i, :) * obj.Risk * PopDec(i, :)';

                PopObj(i, 2) =- sum(PopDec(i, :) * obj.Yield);
            end

            obj.evaluated = obj.evaluated + size(PopDec, 1);

            pause(size(PopObj, 1) * 14.5 + floor(14.5 * 0.05 * rand))
%             varargout = {PopObj};
        end
    end

end
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP7.p
% E:/CODES/DDMOP-master/DDMOP_Exp/Problems\DDMOP7.p
% b'v01.00v00.00\x00\x0e\xa0\x1c\x81\x02_\xb1\x00\x00\x03\xb0\x00\x00\x04A\x00\x00\x07\\'
