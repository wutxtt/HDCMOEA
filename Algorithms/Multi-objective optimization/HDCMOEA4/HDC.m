function Labels = HDC(varargin)
% Hyperdimentional computing classification

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

persistent model

    if nargin == 2 % 训练模型
        
        %% 索引超向量B初始化
        model.data  = [varargin{1};varargin{2}]; % 输入数据分为好解和坏解 Pgood和Pbad
        model.label = [true(1,size(varargin{1},1)),false(1,size(varargin{2},1))]; %Pgood为true Pbad为false
        [N,D] = size(model.data); % 种群的个数N 解的维数D
        h = 10000; % 超维向量的维数
        rho = 0.05; % 超维向量期望的1值个数
        BaseV = sparseVectorsGenerator_sparsity(D,h,rho);  % 生成D个基超向量Base_Vector
        model.BaseV = BaseV;
        model.h = h;

        %% 近似变量值编码方法 / 变量值超向量初始化
        numF = 50;  % 特征值超向量的个数为10
        FeatureV = zeros(numF,h);  % 初始化numF个特征值超向量 分别对应0.05,0.15,0.25,0.35,...,0.95
        FeatureV = single(FeatureV);
        Feature1 = sparseVectorsGenerator_sparsity(1,h,rho);  % 生成第一个特征值超向量Feature_Vector1
        FeatureV(1,:) = Feature1;
        numflip = h/numF;  % 生成余下9个特征值超向量每次变换的位数
        for i = 1:numF-1  % 循环生成余下特征值超向量
            Bit = randi([1,h],1,numflip);   % 随机生成下一个特征值超向量要转变的维度
            newF = FeatureV(i,:); %选用上一个已生成的特征值超向量
            temp = ~newF;
            newF(1,Bit) = temp(1,Bit); % 将所选维度倒置
            FeatureV(i+1,:) = newF;
        end
        
        model.Feature = FeatureV;
        %% 决策变量正则化
%         printf(FeatureV);
        nPgood = size(varargin{1},1);
        Pgoodmax = max(varargin{1});
        Pgoodmin = min(varargin{1});
        NormPgood = (varargin{1}-repmat(Pgoodmin,nPgood,1))./repmat(Pgoodmax-Pgoodmin,nPgood,1);  % 标准化Pgood到[0,1]
        
        nPbad = size(varargin{2},1);
        Pbadmax = max(varargin{2});
        Pbadmin = min(varargin{2});
        NormPbad = (varargin{2}-repmat(Pbadmin,nPbad,1))./repmat(Pbadmax-Pbadmin,nPbad,1);  % 标准化Pbad到[0,1]
        
        %% "非支配解"类别解的编码
        HyperVector = zeros(N,h);
        PFeature = zeros(D,h); % 初始化一个解的特征值超向量
        Value_Feature = linspace(0.05,0.95,numF);  % 创建近似特征值的等差数列
        model.Feature_value = Value_Feature;
        for i = 1:nPgood
            for j = 1:D
                [closestVALUE, closestPOSITION]=min(abs(Value_Feature-NormPgood(i,j)));
%                 PFeature(i,:) = FeatureV(closestPOSITION);
                HyperVariable = xor(BaseV(j,:),FeatureV(closestPOSITION,:));
                HyperVector(i,:) = HyperVector(i,:) + single(HyperVariable(1,:));
            end
        end
        GoodHyperVector = sum(HyperVector);
        model.GH = GoodHyperVector;
        
        %% "支配解"类别解的编码
        HyperVector = zeros(N,h);
        for i = 1:nPbad
            for j = 1:D
                [closestVALUE, closestPOSITION]=min(abs(Value_Feature-NormPbad(i,j)));
%                 PFeature(i,:) = FeatureV(closestPOSITION);
                HyperVariable = xor(BaseV(j,:),FeatureV(closestPOSITION,:));
                HyperVector(i,:) = HyperVector(i,:) + single(HyperVariable(1,:));
            end
        end
        BadHyperVector = sum(HyperVector);
        model.BH = BadHyperVector;
        
        

    else % Predict
        %% 基于HDC分类模型的预选择
        [N,D] = size(varargin{1});
        Pmax = max(varargin{1});
        Pmin = min(varargin{1});
        NormP = (varargin{1}-repmat(Pmin,N,1))./repmat(Pmax-Pmin,N,1);  % 标准化所有候选解到[0,1]
        BaseV = model.BaseV;
        FeatureV = model.Feature;
        Value_Feature = model.Feature_value;
        HyperVector = zeros(N,model.h);
        Labels = zeros(N,1);
        for i = 1:N  % 对每个生成的候选解分配标签
            for j = 1:D
                [closestVALUE, closestPOSITION]=min(abs(Value_Feature-NormP(i,j))); % 确定该维度特征值对应的变量值超向量
%                 PFeature(i,:) = FeatureV(closestPOSITION);
                HyperVariable = xor(BaseV(j,:),FeatureV(closestPOSITION,:)); % 索引超向量与变量值超向量与或
                HyperVector(i,:) = HyperVector(i,:) + single(HyperVariable(1,:)); % 所有维度的超向量按位相加
            end
            cosdisGood = pdist2(HyperVector(i,:),model.GH,'cosine');  % 计算该候选解与"非支配解"类别超向量的余弦相似度
            cosdisBad = pdist2(HyperVector(i,:),model.BH,'cosine');  % 计算该候选解与"支配解"类别超向量的余弦相似度
            if cosdisGood < cosdisBad
                Labels(i) = 1;
            else
                Labels(i) = 0;
            end
            
        end
        
%         Distance = pdist2(varargin{1},model.data);
%         [~,rank] = sort(Distance,2);
%         Labels   = sum(model.label(rank(:,1:5))==1,2) > 2;
    end
end