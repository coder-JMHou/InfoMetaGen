function problem = build_hologram_problem(known_abs_spatial, metrix, varargin)
%BUILD_HOLOGRAM_PROBLEM
% 构造一个可直接喂给 MATLAB 内置优化器的全局优化问题
%
% 用途：
%   problem = build_hologram_problem(known_abs_spatial, metrix, ...)
%
% 返回：
%   problem.objFcn      : 目标函数句柄
%   problem.evaluate    : 给定 x 后，完整前向计算并输出结果
%   problem.decode      : 将优化变量解码为相位/编码
%   problem.lb, ub      : 上下界
%   problem.x0          : 初值
%   problem.intcon      : 整数变量索引（GA整数优化时使用）
%   problem.nvars       : 变量维度（这里就是 num_x*num_y）
%   problem.A, b, Aeq, beq, nonlcon : 空约束，方便直接传给优化器
%
% 支持的设计变量模式 designMode:
%   'integer8'         : 变量直接是 0~7 的整数编码（推荐给 GA）
%   'continuous8quant' : 变量是 [0,2pi] 连续相位，但在目标函数内部量化到8态
%   'continuous'       : 变量是 [0,2pi] 连续相位，不量化
%
% 说明：
%   1) 这是“全 64x64 直接优化”框架，不做任何分块
%   2) 默认目标函数：
%        score = rmse + beta*(1-eta_target) + gamma*(1-pcc_pos)
%      其中：
%        rmse      = 重建图像与目标图像的均方根误差
%        eta_target= 目标区域能量占比
%        pcc_pos   = max(PCC,0)
%
% 作者建议：
%   - 若最终是离散8态超表面，GA + integer8 最符合问题本质
%   - PSO / SA / patternsearch 推荐 continuous8quant
%

%% =========================
% 输入解析
%% =========================
p = inputParser;

addRequired(p, 'known_abs_spatial', @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'metrix', @(x) isnumeric(x) && ismatrix(x));

addParameter(p, 'horn_imamp_single', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'optimizer', 'ga', @(x) ischar(x) || isstring(x));
addParameter(p, 'designMode', '', @(x) ischar(x) || isstring(x));

addParameter(p, 'quantLevel', 8, @(x) isnumeric(x) && isscalar(x) && x>=2);
addParameter(p, 'beta', 0.25, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p, 'gamma', 0.00, @(x) isnumeric(x) && isscalar(x) && x>=0);
addParameter(p, 'maskThreshold', 0.05, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<1);

addParameter(p, 'x0', [], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'rngSeed', 2026, @(x) isnumeric(x) && isscalar(x));

parse(p, known_abs_spatial, metrix, varargin{:});
R = p.Results;

%% =========================
% 基本尺寸
%% =========================
known_abs_spatial = double(known_abs_spatial);
[num_y, num_x] = size(known_abs_spatial);
nvars = num_x * num_y;

if size(metrix,2) ~= nvars
    error('metrix 列数必须等于 num_x*num_y');
end
if size(metrix,1) ~= nvars
    warning('当前默认假设观测面也是 num_x*num_y。若不是，请自行修改 reshape 相关逻辑。');
end

%% =========================
% 入射场
%% =========================
if isempty(R.horn_imamp_single)
    horn_imamp_single = ones(nvars,1);
else
    horn_imamp_single = R.horn_imamp_single(:);
end

if numel(horn_imamp_single) ~= nvars
    error('horn_imamp_single 长度必须等于 num_x*num_y');
end

%% =========================
% 自动选择 designMode
%% =========================
optimizer = lower(string(R.optimizer));

if strlength(string(R.designMode)) == 0
    if optimizer == "ga"
        designMode = "integer8";
    else
        designMode = "continuous8quant";
    end
else
    designMode = lower(string(R.designMode));
end

validModes = ["integer8","continuous8quant","continuous"];
if ~any(designMode == validModes)
    error('designMode 必须是 integer8 / continuous8quant / continuous');
end

%% =========================
% 目标图像归一化与 mask
%% =========================
target_ref_2d = normalize01_local(known_abs_spatial);

mask2d = target_ref_2d > R.maskThreshold * max(target_ref_2d(:));
if ~any(mask2d(:))
    mask2d = target_ref_2d > 0;
end
if ~any(mask2d(:))
    mask2d = true(size(target_ref_2d));
end

%% =========================
% 优化变量边界与初值
%% =========================
% rng(R.rngSeed);

switch designMode
    case "integer8"
        lb = zeros(1, nvars);
        ub = (R.quantLevel - 1) * ones(1, nvars);
        intcon = 1:nvars;

        if isempty(R.x0)
            x0 = randi([0, R.quantLevel-1], 1, nvars);
        else
            x0 = round(R.x0(:)).';
            if numel(x0) ~= nvars
                error('x0 长度必须等于 num_x*num_y');
            end
            x0 = min(max(x0, lb), ub);
        end

    case "continuous8quant"
        lb = zeros(1, nvars);
        ub = 2*pi * ones(1, nvars);
        intcon = [];

        if isempty(R.x0)
            x0 = 2*pi * rand(1, nvars);
        else
            x0 = R.x0(:).';
            if numel(x0) ~= nvars
                error('x0 长度必须等于 num_x*num_y');
            end
            x0 = mod(x0, 2*pi);
        end

    case "continuous"
        lb = zeros(1, nvars);
        ub = 2*pi * ones(1, nvars);
        intcon = [];

        if isempty(R.x0)
            x0 = 2*pi * rand(1, nvars);
        else
            x0 = R.x0(:).';
            if numel(x0) ~= nvars
                error('x0 长度必须等于 num_x*num_y');
            end
            x0 = mod(x0, 2*pi);
        end
end

%% =========================
% 输出 problem 结构体
%% =========================
problem.num_x = num_x;
problem.num_y = num_y;
problem.nvars = nvars;

problem.optimizer = char(optimizer);
problem.designMode = char(designMode);
problem.quantLevel = R.quantLevel;

problem.target_ref_2d = target_ref_2d;
problem.mask2d = mask2d;
problem.horn_imamp_single = horn_imamp_single;
problem.metrix = metrix;

problem.beta = R.beta;
problem.gamma = R.gamma;

problem.lb = lb;
problem.ub = ub;
problem.x0 = x0;
problem.intcon = intcon;

problem.A = [];
problem.b = [];
problem.Aeq = [];
problem.beq = [];
problem.nonlcon = [];

problem.objFcn = @objectiveFcn;
problem.evaluate = @evaluateSolution;
problem.decode = @decodeVariable;

%% =========================
% ========== 嵌套函数 ==========
%% =========================

    function score = objectiveFcn(x)
        evalRes = forwardAndMetrics(x);
        score = evalRes.score;
    end

    function decodeRes = decodeVariable(x)
        x = x(:).';
        if numel(x) ~= nvars
            error('输入变量长度必须等于 num_x*num_y');
        end

        switch designMode
            case "integer8"
                code_vec = round(x);
                code_vec = min(max(code_vec, 0), R.quantLevel-1);

                phase_step = 2*pi / R.quantLevel;
                phase_used_vec = (code_vec + 0.5) * phase_step;
                phase_raw_vec  = phase_used_vec;

            case "continuous8quant"
                phase_raw_vec = mod(x, 2*pi);

                phase_step = 2*pi / R.quantLevel;
                code_vec = floor(phase_raw_vec / phase_step);
                code_vec(code_vec == R.quantLevel) = R.quantLevel - 1;

                phase_used_vec = (code_vec + 0.5) * phase_step;

            case "continuous"
                phase_raw_vec = mod(x, 2*pi);
                phase_used_vec = phase_raw_vec;

                phase_step = 2*pi / R.quantLevel;
                code_vec = floor(phase_raw_vec / phase_step);
                code_vec(code_vec == R.quantLevel) = R.quantLevel - 1;
        end

        decodeRes.x = x;
        decodeRes.phase_raw_vec = phase_raw_vec(:);
        decodeRes.phase_used_vec = phase_used_vec(:);
        decodeRes.code_vec = code_vec(:);

        decodeRes.phase_raw_2d = reshape(phase_raw_vec, num_y, num_x);
        decodeRes.phase_used_2d = reshape(phase_used_vec, num_y, num_x);
        decodeRes.code_2d = reshape(code_vec, num_y, num_x);

        decodeRes.phase_used_deg_2d = rad2deg(mod(decodeRes.phase_used_2d, 2*pi));
        decodeRes.phase_raw_deg_2d  = rad2deg(mod(decodeRes.phase_raw_2d, 2*pi));
    end

    function evalRes = forwardAndMetrics(x)
        dec = decodeVariable(x);

        signal_estimate_spatial = horn_imamp_single .* exp(1i * dec.phase_used_vec);
        target = metrix * signal_estimate_spatial;

        amp_vec = abs(target);
        amp_2d = reshape(amp_vec, num_y, num_x);
        amp_norm_2d = amp_2d ./ (max(amp_2d(:)) + eps);

        ref_2d = target_ref_2d;

        rmse = sqrt(mean((amp_norm_2d(:) - ref_2d(:)).^2));

        power2d = amp_2d.^2;
        eta_target = sum(power2d(mask2d), 'all') / (sum(power2d(:), 'all') + eps);

        power_in = horn_imamp_single.^2;
        eff_total = sum(power2d(:), 'all') / (sum(power_in(:)) + eps);

        pcc = corrcoef_safe(ref_2d(:), amp_norm_2d(:));
        pcc_pos = max(pcc, 0);

        score = rmse + R.beta * (1 - eta_target) + R.gamma * (1 - pcc_pos);

        evalRes.score = score;
        evalRes.rmse = rmse;
        evalRes.eta_target = eta_target;
        evalRes.eff_total = eff_total;
        evalRes.pcc = pcc;
        evalRes.pcc_pos = pcc_pos;

        evalRes.target_complex_vec = target;
        evalRes.target_amp_vec = amp_vec;
        evalRes.target_amp_2d = amp_2d;
        evalRes.target_amp_norm_2d = amp_norm_2d;

        evalRes.decode = dec;
    end

    function result = evaluateSolution(x)
        evalRes = forwardAndMetrics(x);
        dec = evalRes.decode;

        result = struct();
        result.xbest = x(:).';

        result.score = evalRes.score;
        result.rmse = evalRes.rmse;
        result.eta_target = evalRes.eta_target;
        result.eff_total = evalRes.eff_total;
        result.pcc = evalRes.pcc;

        result.phase_raw_2d = dec.phase_raw_2d;
        result.phase_used_2d = dec.phase_used_2d;
        result.phase_raw_deg_2d = dec.phase_raw_deg_2d;
        result.phase_used_deg_2d = dec.phase_used_deg_2d;
        result.code_2d = dec.code_2d;

        result.target_amp_2d = evalRes.target_amp_2d;
        result.target_amp_norm_2d = evalRes.target_amp_norm_2d;
        result.target_complex_vec = evalRes.target_complex_vec;

        result.target_ref_2d = target_ref_2d;
        result.mask2d = mask2d;
    end

end

%% =========================
% 局部函数
%% =========================
function X = normalize01_local(X)
X = double(X);
xmin = min(X(:));
xmax = max(X(:));
if abs(xmax - xmin) < eps
    X = zeros(size(X));
else
    X = (X - xmin) / (xmax - xmin);
end
end

function pcc = corrcoef_safe(a, b)
a = a(:);
b = b(:);

if std(a) < eps || std(b) < eps
    pcc = 0;
    return;
end

R = corrcoef(a, b);
if numel(R) < 4 || any(isnan(R(:)))
    pcc = 0;
else
    pcc = R(1,2);
end
end