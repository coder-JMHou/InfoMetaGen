function result = optimize_hologram_builtin(problem, optimizer, solverOptions)
%OPTIMIZE_HOLOGRAM_BUILTIN
% 用统一接口调用 MATLAB 自带优化器
%
% 用法：
%   result = optimize_hologram_builtin(problem, 'ga', opts)
%   result = optimize_hologram_builtin(problem, 'particleswarm', opts)
%   result = optimize_hologram_builtin(problem, 'simulannealbnd', opts)
%   result = optimize_hologram_builtin(problem, 'patternsearch', opts)
%
% 输入：
%   problem       : build_hologram_problem 返回的结构体
%   optimizer     : 'ga' / 'particleswarm' / 'simulannealbnd' / 'patternsearch'
%   solverOptions : 对应优化器的 options；若为空则自动给默认值
%
% 输出：
%   result 中包含：
%       xbest, fbest, exitflag, output
%       以及 problem.evaluate(xbest) 的所有结果
%

if nargin < 2 || isempty(optimizer)
    optimizer = problem.optimizer;
end
if nargin < 3
    solverOptions = [];
end

optimizer = lower(string(optimizer));

switch optimizer
    case "ga"
        if isempty(solverOptions)
            solverOptions = optimoptions('ga', ...
                'Display', 'iter', ...
                'PopulationSize', 30, ...
                'MaxGenerations', 20, ...
                'UseParallel', false);
        end

        if isempty(problem.intcon)
            [xbest, fbest, exitflag, output, population, scores] = ga( ...
                problem.objFcn, problem.nvars, ...
                problem.A, problem.b, problem.Aeq, problem.beq, ...
                problem.lb, problem.ub, problem.nonlcon, ...
                solverOptions);
        else
            [xbest, fbest, exitflag, output, population, scores] = ga( ...
                problem.objFcn, problem.nvars, ...
                problem.A, problem.b, problem.Aeq, problem.beq, ...
                problem.lb, problem.ub, problem.nonlcon, ...
                problem.intcon, solverOptions);
        end

        evalRes = problem.evaluate(xbest);

        result = evalRes;
        result.optimizer = char(optimizer);
        result.fbest = fbest;
        result.exitflag = exitflag;
        result.output = output;
        result.population = population;
        result.scores = scores;

    case "particleswarm"
        if ~isempty(problem.intcon)
            error(['particleswarm 不能直接做整数约束。', ...
                   '请把 designMode 改成 ''continuous8quant'' 或 ''continuous''。']);
        end

        if isempty(solverOptions)
            solverOptions = optimoptions('particleswarm', ...
                'Display', 'iter', ...
                'SwarmSize', 20, ...
                'MaxIterations', 20, ...
                'UseParallel', false);
        end

        [xbest, fbest, exitflag, output] = particleswarm( ...
            problem.objFcn, problem.nvars, problem.lb, problem.ub, solverOptions);

        evalRes = problem.evaluate(xbest);

        result = evalRes;
        result.optimizer = char(optimizer);
        result.fbest = fbest;
        result.exitflag = exitflag;
        result.output = output;

    case "simulannealbnd"
        if ~isempty(problem.intcon)
            error(['simulannealbnd 不能直接做整数约束。', ...
                   '请把 designMode 改成 ''continuous8quant'' 或 ''continuous''。']);
        end

        if isempty(solverOptions)
            solverOptions = optimoptions('simulannealbnd', ...
                'Display', 'iter', ...
                'MaxIterations', 200, ...
                'FunctionTolerance', 1e-6);
        end

        [xbest, fbest, exitflag, output] = simulannealbnd( ...
            problem.objFcn, problem.x0, problem.lb, problem.ub, solverOptions);

        evalRes = problem.evaluate(xbest);

        result = evalRes;
        result.optimizer = char(optimizer);
        result.fbest = fbest;
        result.exitflag = exitflag;
        result.output = output;

    case "patternsearch"
        if ~isempty(problem.intcon)
            error(['patternsearch 通常不直接做整数约束。', ...
                   '请把 designMode 改成 ''continuous8quant'' 或 ''continuous''。']);
        end

        if isempty(solverOptions)
            solverOptions = optimoptions('patternsearch', ...
                'Display', 'iter', ...
                'UseCompletePoll', true, ...
                'UseCompleteSearch', false);
        end

        [xbest, fbest, exitflag, output] = patternsearch( ...
            problem.objFcn, problem.x0, ...
            problem.A, problem.b, problem.Aeq, problem.beq, ...
            problem.lb, problem.ub, problem.nonlcon, ...
            solverOptions);

        evalRes = problem.evaluate(xbest);

        result = evalRes;
        result.optimizer = char(optimizer);
        result.fbest = fbest;
        result.exitflag = exitflag;
        result.output = output;

    otherwise
        error('暂不支持该优化器：%s', optimizer);
end
end