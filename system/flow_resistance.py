from abc import ABC
import numpy as np
import system.interpolation as ip


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

    def update(self):
        pass


class ConstantFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        self.value = zeta_dict['value']


class WallFrictionFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        self.channel = channel
        self.method = zeta_dict.get('method', 'Blasius')
        self.type = zeta_dict.get('type', 'Darcy')
        self.value = np.zeros(self.channel.n_ele)

    def update(self):
        f = 1.0
        if type == 'Darcy':
            f = 4.0
        elif type == 'Fanning':
            f = 1.0
        else:
            ValueError('friction factor type can only be Darcy or Fanning')
        lam = np.zeros(self.value.shape)
        turb = np.zeros(self.value.shape)
        reynolds = ip.interpolate_1d(self.channel.reynolds)
        if self.method == 'Blasius':
            lam = np.divide(f * 16.0, reynolds, out=lam, where=reynolds > 0.0)
            turb = f * 0.0791 * np.power(reynolds, -0.25, out=turb,
                                         where=reynolds > 0.0)
            factor = np.where(reynolds < 2200.0, lam, turb)
            # factor = lam
            self.value[:] = self.channel.dx / self.channel.d_h * factor
        else:
            raise NotImplementedError


class JunctionFlowResistance(FlowResistance):
    def __init__(self, channel, zeta_dict, **kwargs):
        self.channel = channel
        self.zeta_const = zeta_dict.get('value', 0.0)
        self.factor = zeta_dict['factor']
        self.value = np.zeros(self.channel.n_ele)

    def update(self):
        ref_velocity = np.max(self.channel.velocity)
        self.value[:] = self.zeta_const
        if np.abs(ref_velocity > 0.0):
            self.value[:] += \
                self.factor * np.log(self.channel.velocity[:-1] / ref_velocity)


