classdef HDCMOEA5 < ALGORITHM
% <multi> <real> <expensive>
% Classification and Pareto domination based multi-objective evolutionary
% algorithm
% M --- 3 --- Number of generated offsprings for each solution

%------------------------------- Reference --------------------------------
% J. Zhang, A. Zhou, and G. Zhang, A classification and Pareto domination
% based multiobjective evolutionary algorithm, Proceedings of the IEEE
% Congress on Evolutionary Computation, 2015, 2883-2890.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            M = Algorithm.ParameterSet(3);

            %% Generate random population
            Population   = Problem.Initialization(); % 种群初始化
            A = Population; % 存档集初始化
            [Pgood,Pbad] = NDS(Population,floor(Problem.N/2)); %存档集划分

            %% Optimization
            while Algorithm.NotTerminated(Population)
                HDC(Pgood.decs, Pbad.decs); % 分类器HDC构建
%                 KNN(Pgood.decs,Pbad.decs);
                Offspring  = Operator(Problem,Population,M); % 候选解生成 && 基于HDC分类模型的预选择 && 评价好的新子代解
                A = [A, Offspring]; % 存档集更新
                [Pgood,Pbad] = NDS(A,floor(Problem.N/2)); %存档集划分
                Population = NDS(A,Problem.N);
                if Problem.FE >= Problem.maxFE % 停止准则 从存档集中输出非支配解
%                     [FrontNo,~] = NDSort(A.objs,A.cons,Problem.N);
%                     Algorithm.NotTerminated(A(FrontNo==1));
                    P = NDS(A,Problem.N);
                    Algorithm.NotTerminated(P);
                end 
            end
        end
    end
end