"""连续筛查使用的分布函数（仅依赖标准库）。"""
import math

_EPS = 1e-15


def _betacf(a, b, x):
    """连分式展开的补不全 Beta 函数，Numerical Recipes 口径。"""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < _EPS:
        d = _EPS
    d = 1.0 / d
    result = d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < _EPS:
            d = _EPS
        c = 1.0 + aa / c
        if abs(c) < _EPS:
            c = _EPS
        d = 1.0 / d
        result *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < _EPS:
            d = _EPS
        c = 1.0 + aa / c
        if abs(c) < _EPS:
            c = _EPS
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < 3e-14:
            break
    return result


def incomplete_beta_regularized(x, a, b):
    """正则化不全 Beta 函数 I_x(a, b)，x 必须位于 [0, 1]。"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_prefix = (a * math.log(x) + b * math.log1p(-x)
                  + math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b))
    prefix = math.exp(log_prefix)
    if x < (a + 1.0) / (a + b + 2.0):
        return prefix * _betacf(a, b, x) / a
    return 1.0 - prefix * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(statistic, df):
    """Student t 统计量的双侧 p 值；df 为正整数。"""
    if df <= 0:
        raise ValueError('自由度必须为正整数')
    x = df / (df + statistic * statistic)
    return incomplete_beta_regularized(x, df / 2.0, 0.5)


def lag1_autocorrelation(values):
    """滞后 1 阶 Pearson 自相关；样本不足时返回 None。"""
    if len(values) < 3:
        return None
    mean = math.fsum(values) / len(values)
    centered = [value - mean for value in values]
    denominator = math.fsum(z * z for z in centered)
    if denominator <= 0.0:
        return 0.0
    numerator = math.fsum(centered[i] * centered[i - 1]
                          for i in range(1, len(centered)))
    return numerator / denominator
