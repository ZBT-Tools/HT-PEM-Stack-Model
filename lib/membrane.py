import numpy as np
from abc import ABC, abstractmethod
from . import layers as layers
from data import global_parameters as g_par


class Membrane(ABC, layers.SolidLayer):
    def __new__(cls, membrane_dict, dx, **kwargs):
        model_type = membrane_dict.get('type', 'Constant')
        if model_type is 'Constant':
            return super(Membrane, cls).__new__(Constant)
        elif model_type is 'Springer':
            return super(Membrane, cls).__new__(SpringerMembrane)
        elif model_type is 'Kvesic':
            return super(Membrane, cls).__new__(KvesicMembrane)
        else:
            raise NotImplementedError('Specified membrane model not '
                                      'implemented. Available models are '
                                      'Constant, Springer, and Kvesic.')

    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx)

        self.temp = np.zeros(self.dx.shape)
        # membrane temperature
        self.ionic_conductivity = \
            membrane_dict.get('ionic conductivity', 1.0e-3)
        # constant ionic conductivity of membrane
        self.omega_ca = np.zeros(self.dx.shape)
        # area specific membrane resistance
        self.omega = np.zeros(self.dx.shape)
        # membrane resistance
        self.calc_loss = membrane_dict.get('calc_loss', True)

        self.v_loss = np.zeros(self.dx.shape)
        # voltage loss at the membrane
        self.ionic_conductance = self.calc_conductance(self.ionic_conductivity)

    @abstractmethod
    def calc_ionic_resistance(self, *args):
        pass

    def calc_voltage_loss(self, current_density):
        """
        Calculates the voltage loss at the membrane.
        """
        if not self.calc_loss:
            self.v_loss[:] = 0.
        else:
            self.v_loss[:] = self.omega_ca * current_density

    def update(self, current_density, humidity, *args):
        self.calc_ionic_resistance(humidity, *args)
        self.calc_voltage_loss(current_density)


class WaterTransportMembrane(Membrane, ABC):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)

        self.vapour_coeff = membrane_dict['vapour transport coefficient']
        self.acid_group_conc = membrane_dict['acid group concentration']
        self.w_cross_flow = np.zeros(self.dx.shape)
        self.faraday_const = g_par.constants['F']
        # water cross flux through the membrane

    def calc_cross_water_flux(self, current_density, humidity):
        """
        Calculates the water cross flux through the membrane
        according to (Springer, 1991).
        """
        dw = 2.1e-7 * np.exp(-2436. / self.temp)
        water_content = 0.043 + 17.81 * humidity \
            - 39.85 * humidity ** 2. + 36. * humidity ** 3.

        divisor = 2. * self.vapour_coeff * self.acid_group_conc \
            * self.faraday_const
        zeta_plus = \
            water_content[0] + water_content[1] + current_density / divisor
        zeta_negative = (water_content[0] - water_content[1]
                         + 5. * current_density / divisor) \
            / (1. + dw * zeta_plus / (self.thickness * self.vapour_coeff))
        m_c = 0.5 * (zeta_plus + zeta_negative)
        m_a = 0.5 * (zeta_plus - zeta_negative)
        self.w_cross_flow[:] = \
            current_density / self.faraday_const + self.acid_group_conc \
            * dw * (m_a ** 2. - m_c ** 2.) / (2. * self.thickness)

    def update(self, current_density, humidity, *args):
        self.calc_cross_water_flux(current_density, humidity)
        super().update(current_density, humidity, *args)


class Constant(Membrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        self.w_cross_flow = np.zeros_like(self.dx)
        # water cross flux through the membrane
        self.omega[:] = 1.0 / self.ionic_conductance[0]
        self.omega_ca[:] = self.omega * self.area_dx

    def calc_ionic_resistance(self, *args):
        pass


class SpringerMembrane(WaterTransportMembrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)

    def calc_ionic_resistance(self, humidity):
        """
        Calculates the membrane resistivity
        for NT-PEMFC according to (Springer, 1991).
        """
        humidity = np.average(humidity, axis=0)
        lambda_springer = \
            np.where(humidity < 1.0,
                     0.043 + 17.81 * humidity
                     - 39.85 * humidity ** 2. + 36. * humidity ** 3.,
                     14.0 + 1.4 * (humidity - 1.0))
        lambda_springer[lambda_springer < 1.0] = 1.0
        mem_cond = (0.005139 * lambda_springer - 0.00326) \
            * np.exp(1268.0 * (0.0033 - 1. / self.temp)) * 1e2
        self.omega_ca[:] = self.thickness / mem_cond  # * 1.e-4
        self.omega[:] = self.omega_ca / self.area_dx
        return self.omega, self.omega_ca


class KvesicMembrane(Membrane):
    def __init__(self, membrane_dict, dx, **kwargs):
        super().__init__(membrane_dict, dx, **kwargs)
        self.basic_resistance = membrane_dict['basic resistance']
        # basic electrical resistance of the membrane
        self.temp_coeff = membrane_dict['temperature coefficient']
        # thermal related electrical resistance gain of the membrane

    def calc_ionic_resistance(self, *args):
        self.omega_ca[:] = \
            (self.basic_resistance - self.temp_coeff * self.temp)  # * 1e-2
        self.omega[:] = self.omega_ca / self.area_dx
        return self.omega, self.omega_ca
