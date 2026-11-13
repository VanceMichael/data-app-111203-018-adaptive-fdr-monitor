"""在线多重检验错误预算规则。

本模块实现两种具有明确适用假设的顺序规则：

* ``lord3``：Ramdas 等 (2017, 2018) 与 onlineFDR 包口径的 LORD3。
  适用假设：各次检验的 p 值相互独立（更一般地满足 PRDS 正向依赖），
  且每个真零假设下 p 值为超均匀（super-uniform）分布。在此假设下
  在线 FDR 不超过包络 alpha。使用可求和序列 gamma_t = C / t**2，
  C = 6/pi**2，初始财富 W0 = alpha / 2（满足 W0 <= alpha 且
  sum gamma_t <= 1 的定理条件）。
* ``online_bonferroni``：阈值 beta_t = alpha * gamma_t 的在线 Bonferroni
  包络。除“真零假设下 p 值超均匀”外不要求任何依赖结构；按并集界
  （union bound）控制 FWER <= alpha，因此也控制 FDR <= alpha。
  当序列依赖无法确认满足 LORD3 条件时使用此保守降级规则。

两条规则都是纯顺序的：阈值只取决于当前位置与历史拒绝位置，后来的
记录不可能改变已经花费预算的历史位置。
"""
from dataclasses import dataclass, field
from enum import Enum
from math import pi

GAMMA_C = 6.0 / (pi ** 2)
LORD_INITIAL_FRACTION = 0.5
# 财富下限：理论上财富只趋近于 0，但浮点累减会下溢到精确 0，
# 导致长 horizon 下后续阈值永久归零。钳到此正值保持规则的序贯语义。
WEALTH_FLOOR = 1e-12


def gamma_at(position):
    """可求和权重序列 gamma_t = (6/pi^2) / t^2，t >= 1。"""
    if position < 1:
        raise ValueError('位置必须从 1 开始')
    return GAMMA_C / (position * position)


class Stream(str, Enum):
    LORD3 = 'lord3'
    ONLINE_BONFERRONI = 'online_bonferroni'


@dataclass
class ControllerStatus:
    stream: str
    envelope_alpha: float
    position: int
    wealth: float
    reject_count: int
    reject_positions: list

    def to_dict(self):
        return {'stream': self.stream, 'envelope_alpha': self.envelope_alpha,
                'position': self.position, 'wealth': self.wealth,
                'reject_count': self.reject_count,
                'reject_positions': list(self.reject_positions)}

    @classmethod
    def from_dict(cls, payload):
        return cls(stream=payload['stream'],
                   envelope_alpha=payload['envelope_alpha'],
                   position=payload['position'], wealth=payload['wealth'],
                   reject_count=payload['reject_count'],
                   reject_positions=list(payload['reject_positions']))


@dataclass
class OnlineFDRController:
    envelope_alpha: float
    stream: Stream
    position: int = 0
    wealth: float = field(init=False)
    reject_positions: list = field(default_factory=list)

    def __post_init__(self):
        if not 0.0 < self.envelope_alpha <= 1.0:
            raise ValueError('包络 alpha 必须位于 (0, 1]')
        self.stream = Stream(self.stream)
        self.wealth = LORD_INITIAL_FRACTION * self.envelope_alpha

    @property
    def initial_wealth(self):
        return LORD_INITIAL_FRACTION * self.envelope_alpha

    @property
    def reject_count(self):
        return len(self.reject_positions)

    def threshold(self):
        """返回下一次检验（position + 1）的预算阈值。"""
        return self._threshold_at(self.position + 1)

    def submit(self, p_value):
        """顺序提交一个 p 值，返回本次预算流水并推进位置。"""
        if not 0.0 <= p_value <= 1.0:
            raise ValueError('p 值必须位于 [0, 1]')
        self.position += 1
        current_position = self.position
        alpha = self._threshold_at(current_position)
        rejected = p_value <= alpha
        entry = {'position': current_position, 'stream': self.stream.value,
                 'envelope_alpha': self.envelope_alpha,
                 'alpha': alpha, 'p_value': p_value,
                 'rejected': rejected, 'wealth_before': self.wealth}
        if self.stream is Stream.LORD3:
            if rejected:
                self.wealth += self.envelope_alpha - alpha
                self.reject_positions.append(current_position)
            else:
                self.wealth -= alpha
            self.wealth = max(WEALTH_FLOOR, self.wealth)
            entry['wealth_after'] = self.wealth
        else:
            entry['wealth_after'] = None
        return entry

    def _threshold_at(self, current_position):
        if self.stream is Stream.ONLINE_BONFERRONI:
            return min(self.envelope_alpha,
                       self.envelope_alpha * gamma_at(current_position))
        if current_position == 1:
            return min(self.envelope_alpha, self.wealth)
        alpha_t = self.wealth * gamma_at(current_position)
        for reject_position in self.reject_positions:
            alpha_t += (self.envelope_alpha - self.initial_wealth) * (
                gamma_at(current_position - reject_position))
        return min(self.envelope_alpha, max(0.0, alpha_t))

    def status(self):
        return ControllerStatus(stream=self.stream.value,
                                envelope_alpha=self.envelope_alpha,
                                position=self.position, wealth=self.wealth,
                                reject_count=self.reject_count,
                                reject_positions=list(self.reject_positions))

    @classmethod
    def restore(cls, snapshot):
        status = ControllerStatus.from_dict(snapshot)
        controller = cls(envelope_alpha=status.envelope_alpha,
                         stream=Stream(status.stream))
        controller.position = status.position
        controller.wealth = max(WEALTH_FLOOR, status.wealth)
        controller.reject_positions = list(status.reject_positions)
        return controller
