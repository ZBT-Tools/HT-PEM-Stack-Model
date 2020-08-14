# general imports
from abc import ABC
import numpy as np

# local module imports
from . import interpolation as ip


class FlowResistance(ABC):
    def __new__(cls, channel, zeta_dict, **kwargs):
        zeta_type = zeta_dict.get('type', 'Constant')
        if zeta_type == 'Constant':
            return super(FlowResistance, cls).\
                __new__(ConstantFlowResistance)
        elif zeta_type == 'WallFriction':
            return super(FlowResistance, cls).\
                __new__(WallFrictionFlowResistance)
        elif zeta_type == 'Junction':
            return super(FlowResistance, cls).\
                __new__(JunctionFlowResistance)
        else:
            raise NotImplementedError

    def __init__(self, channel, zeta_dict, **kwargs):
        self.channel = channel
        self.value = None

    def update(self):
        pass


class ConstantFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.value = zeta_dict['value']


class WallFrictionFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.method = zeta_dict.get('method', 'Blasius')
        self.value = np.zeros(self.channel.n_nodes)

    def update(self):
        # reynolds = ip.interpolate_1d(self.channel.reynolds)
        reynolds = self.channel.reynolds
        lam = reynolds < 2200.0
        turb = np.invert(lam)
        lam_id = np.where(lam)
        turb_id = np.where(turb)

        # Darcy friction factor
        factor = np.zeros(reynolds.shape)
        if np.any(lam):
            reynolds_lam = reynolds[lam_id]
            if self.channel.aspect_ratio == 1.0:
                f_reynolds = 64.0
            elif self.channel.cross_shape == 'rectangular':
                eps = self.channel.aspect_ratio
                f_reynolds = 4.0 * 24.0 \
                    / ((1.0 + eps) ** 2.0
                       * (1.0 - (192.0 * eps / np.pi ** 5.0
                                 * np.tanh(np.pi / (2.0 * eps)))))
            else:
                raise NotImplementedError

            factor[lam_id] = \
                np.divide(f_reynolds, reynolds_lam, where=reynolds_lam > 0.0)

        if np.any(turb):
            reynolds_turb = reynolds[turb_id]
            if self.method == 'Blasius':
                factor[turb_id] = 0.3164 * np.power(reynolds_turb, -0.25)
            else:
                raise NotImplementedError
        np.seterr(under='ignore')
        self.value[:] = self.channel.dx_node / self.channel.d_h * factor
        np.seterr(under='raise')


class JunctionFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        super().__init__(channel, zeta_dict, **kwargs)
        self.zeta_const = zeta_dict.get('value', 0.0)
        self.factor = zeta_dict['factor']
        self.value = np.zeros(self.channel.n_nodes)

    def update(self):
        ref_velocity = np.max(self.channel.velocity)
        self.value[:] = self.zeta_const
        if np.abs(ref_velocity > 0.0):
            self.value[:] += \
                self.factor * np.log(self.channel.velocity / ref_velocity)
