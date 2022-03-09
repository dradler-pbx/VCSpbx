import numpy as np
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI as CPPSI
import CoolProp.CoolProp as CP


def LMTD_calc(Thi, Tho, Tci, Tco):
    # calculating the logaritmic mean temperature of two fluids with defined "hot" and "cold" fluid
    dT1 = Thi - Tco
    dT2 = Tho - Tci
    if dT2 == 0:
        dT2 = 0.01
    if dT1 == dT2:
        LMTD = dT1
    else:
        LMTD = (dT1 - dT2) / np.log(dT1 / dT2)
    if dT1 < 0:
        return 0.0
    # prevent NaN values:
    if np.isnan(LMTD):
        LMTD = 1e-6
    return LMTD

def dhCOND(TC, medium):
    Tcrit = CPPSI('TCRIT', medium)
    if TC > Tcrit:
        dh = -(TC - Tcrit) * 1000
    else:
        dh = CPPSI("H", "T", TC, "Q", 1, medium) - CPPSI("H", "T", TC, "Q", 0, medium)
    return dh


class System:
    def __init__(self, id: str, tolerance: float, n_max: int = 100, fun_tol: float = 0.1):
        self.id = id
        self.components = []
        self.junctions = []
        self.params = {}
        self.tolerance = tolerance

        self.n_max = n_max
        self.fun_tol = fun_tol

        self.residual_enthalpy = None
        self.residual_functions = {}

    def run(self, full_output=False):
        # initialize the system and all components
        self.initialize()
        # first get the current enthalpy values
        old_enthalpies = self.get_junction_enthalpies()
        counter = 0
        while True:
            for comp in self.components:
                comp.calc()
            # get the new enthalpy values
            new_enthalpies = self.get_junction_enthalpies()

            # calculate the delta
            abs_delta = np.abs(old_enthalpies - new_enthalpies)

            # get the residual of the functions
            fun_residual = np.array([])
            for comp in self.components:
                res = np.array(comp.get_function_residual())
                fun_residual = np.append(fun_residual, np.abs(res))

            # check if delta is lower than tolerance
            if np.max(abs_delta) < self.tolerance and np.max(fun_residual) < self.fun_tol:
                break
            counter += 1
            if counter > self.n_max:
                print('Reached {} iterations without solution'.format(counter))
                for comp in self.components:
                    self.residual_functions[comp.id] = comp.get_function_residual()
                self.residual_enthalpy = abs_delta
                return False

            old_enthalpies = new_enthalpies.copy()
        # Helper
            for junc in self.junctions:
                print([junc.id, junc.T, junc.h])
            print('--')

        self.residual_enthalpy = abs_delta

        for comp in self.components:
            self.residual_functions[comp.id] = comp.get_function_residual()

        if full_output:
            print('---')
            print('Iteration finished after {} iterations'.format(counter))
            print('Residual enthalpies difference: {}'.format(abs_delta))
            print('---\nJUNCTIONS:')
            for junc in self.junctions:
                print(junc.id)
                print('Pressure: {:2f}bar'.format(junc.p/1e5))
                print('Temperature: {:2f}°C'.format(junc.T-273.15))
                print('Enthalpy: {:2f}J/kg'.format(junc.h))
                print('Massflow: {:2f}kg/s'.format(junc.mdot))
                print('---')

        return True

    def initialize(self):
        for comp in self.components:
            comp.check_junctions()
            comp.initialize()

    def update_parameters(self):
        pass

    def add_component(self, comp):
        self.components.append(comp)

    def add_junction(self, junc):
        self.junctions.append(junc)

    def get_junction_enthalpies(self):
        return np.array([junc.get_enthalpy() for junc in self.junctions])


class Component:
    def __init__(self, id: str, system: object):
        self.system = system
        system.add_component(self)

        self.id = id

        self.parameters = {}

        self.junctions = {'inlet_A': None, 'outlet_A': None}

    def print_parameters(self):
        print(self.parameters)

    def calc(self):
        print('{} has no defined calc() function!'.format(self.id))

    def check_junctions(self):
        for key in self.junctions:
            if not self.junctions[key]:
                raise ValueError('{} has no junction at port {}.'.format(self.id, key))

    def get_function_residual(self):
        return 0.0

class Junction:
    def __init__(self, id: str, system: object, medium: str, upstream_component:object, upstream_id: str, downstream_component: object, downstream_id: str, mdot_init: float, p_init: float, h_init: float):
        self.medium = medium

        self.id = id
        self.system = system
        system.add_junction(self)

        self.set_values(mdot_init, p_init, h_init)

        if upstream_component.junctions[upstream_id]:
            print('{} of component {} overwritten!'.format(upstream_id, upstream_component.id))
        upstream_component.junctions[upstream_id] = self
        if downstream_component.junctions[downstream_id]:
            print('{} of component {} overwritten!'.format(downstream_id, downstream_component.id))
        downstream_component.junctions[downstream_id] = self

    def set_values(self, mdot: float = None, p:float = None, h: float = None):
        if mdot:
            self.mdot = mdot
        if p:
            self.p = p
        if h:
            self.h = h

        self.T = CPPSI('T', 'P', self.p, 'H', self.h, self.medium)

        self.s = CPPSI('S', 'P', self.p, 'H', self.h, self.medium)

        try:
            self.x = self.calculate_x()
        except:
            self.x = None

    def get_pressure(self):
        return self.p

    def get_temperature(self):
        return self.T

    def get_massflow(self):
        return self.mdot

    def get_enthalpy(self):
        return self.h

    def get_entropy(self):
        return self.s

    def get_quality(self):
        return self.x

    def calculate_x(self):
        # this defines h for more than just the two phase region
        h_l = CPPSI('H', 'P', self.p, 'Q', 0, self.medium)
        h_v = CPPSI('H', 'P', self.p, 'Q', 1, self.medium)
        return (self.h - h_l)/(h_v-h_l)

class Compressor_efficiency(Component):
    def __init__(self, id: str, system: object, etaS: float, etaV:float, stroke: float, speed: float):
        super().__init__(id, system)
        self.etaS = etaS
        self.etaV = etaV
        self.stroke = stroke
        self.speed = speed

        self.parameters = {'id': id, 'etaS': etaS, 'etaV': etaV, 'speed': speed}

        self.Tin = np.nan
        self.pin = np.nan
        self.pout = np.nan

    def initialize(self):
        pass

    def calc(self):
        self.Tin = self.junctions['inlet_A'].get_temperature()
        self.pin = self.junctions['inlet_A'].get_pressure()
        self.pout = self.junctions['outlet_A'].get_pressure()

        # compressor model based on efficiency parameters: etaS...isentropic / etaV...volumetric
        rho = CPPSI("D", "P", self.pin, "T", self.Tin, "R290")
        mdot = self.speed / 60 * self.stroke * self.etaV * rho  # mass flow
        hin = CPPSI("H", "T", self.Tin, "P", self.pin, "R290")  # inlet enthalpy
        sin = CPPSI("S", "T", self.Tin, "P", self.pin, "R290")  # inlet entropy
        houtS = CPPSI("H", "S", sin, "P", self.pout, "R290")  # enthalpy at outlet under isentropic conditions
        Pel = mdot * (houtS - hin) / self.etaS  # power input
        hout = Pel / mdot + hin  # real outlet enthalpy
        Tout = CPPSI("T", "P", self.pout, "H", hout, "R290")  # outlet temperature

        self.junctions['outlet_A'].set_values(mdot=mdot, h=hout)


class Condenser(Component):
    def __init__(self, id: str, system: object, k: iter, area: float, subcooling: float, T_air_in: float, mdot_air_in: float):
        super().__init__(id, system)
        if len(k) != 3:
            raise ValueError('k needs to be of length 3, but len(k) = {}'.format(len(k)))
        else:
            self.k = k
        self.area = area
        self.dTSC = subcooling
        # self.parameters = {'UA': self.UA, 'subcooling': self.dTSC}

        self.T_air_in = T_air_in
        self.mdot_air = mdot_air_in

        self.TC = None
        self.TAo_desuperheat = None
        self.TAo_condenser = None
        self.TAo_subcool = None
        self.areafraction_desuperheat = None
        self.areafraction_condenser = None
        self.areafraction_subcool = None
        self.p = None

    def initialize(self):
        self.medium = self.junctions['inlet_A'].medium
        self.p = self.junctions['inlet_A'].get_pressure()
        self.TC = CPPSI('T', 'P', self.p, 'Q', 0, self.medium)
        Tmean = (self.T_air_in + self.TC)/2
        self.TAo_desuperheat = Tmean
        self.TAo_condenser = Tmean
        self.TAo_subcool = Tmean
        self.areafraction_desuperheat = 0.1
        self.areafraction_condenser = 0.8
        self.areafraction_subcool = 0.1

    def model(self, x):
        mR = self.junctions['inlet_A'].get_massflow()
        TRin = self.junctions['inlet_A'].get_temperature()
        mA = self.mdot_air
        TAi = self.T_air_in

        # Boundary for fsolve calculation to not cause the logaritmic mean temperature to generate NaN values (neg. logarithm):
        # The outlet air temperature of the superheat section must not exceed the refrigerant inlet temperature.
        if x[1] > TRin:
            x[1] = TRin - 1e-4
        if x[0] - self.dTSC < TAi:
            x[0] = self.T_air_in + self.dTSC + 1e-4

        # calculate material parameters
        cpR = np.zeros(2)
        cpR[0] = CPPSI("C", "T", x[0], "Q", 1, self.medium)
        cpR[1] = CPPSI("C", "T", x[0], "Q", 0, self.medium)
        cpA = CPPSI("C", "T", self.T_air_in, "P", 1.0e5, "AIR")

        # Calculate the mean logaritmic temperature value for all three sections of the condenser
        LMTD = np.zeros(3)
        LMTD[0] = LMTD_calc(TRin, x[0], self.T_air_in, x[1])
        LMTD[1] = LMTD_calc(x[0], x[0], self.T_air_in, x[2])
        LMTD[2] = LMTD_calc(x[0], x[0] - self.dTSC, self.T_air_in, x[3])

        # Formulation of the equation system as according to fsolve documentation ( 0 = ... ).
        # The equation set  and model definition is documented in the model description.
        dh = dhCOND(x[0], self.medium)
        f = np.zeros(7)
        f[0] = mR * cpR[0] * (TRin - x[0]) - mA * cpA * x[4] * (x[1] - TAi)
        f[1] = mR * cpR[0] * (TRin - x[0]) - self.k[0] * x[4] * self.area * LMTD[0]
        f[2] = mR * dh - self.k[1] * x[5] * self.area * LMTD[1]
        f[3] = mR * dh - mA * cpA * x[5] * (x[2] - TAi)
        f[4] = mR * cpR[1] * self.dTSC - mA * cpA * x[6] * (x[3] - TAi)
        f[5] = mR * cpR[1] * self.dTSC - self.k[2] * x[6] * LMTD[2]
        f[6] = 1 - x[4] - x[5] - x[6]

        return f

    def calc(self):

        x = np.zeros(7)
        x[0] = self.TC
        x[1] = self.TAo_desuperheat
        x[2] = self.TAo_condenser
        x[3] = self.TAo_subcool
        x[4] = self.areafraction_desuperheat
        x[5] = self.areafraction_condenser
        x[6] = self.areafraction_subcool

        x = fsolve(self.model, x0=x, xtol=self.system.fun_tol)

        self.TC = x[0]
        self.TAo_desuperheat = x[1]
        self.TAo_condenser = x[2]
        self.TAo_subcool = x[3]
        self.areafraction_desuperheat = x[4]
        self.areafraction_condenser = x[5]
        self.areafraction_subcool = x[6]

        self.p = CPPSI('P', 'T', self.TC, 'Q', 0, self.medium)

        if self.dTSC == 0:
            hout = CPPSI('H', 'P', self.p, 'Q', 0, self.medium)
        else:
            hout = CPPSI('H', 'P', self.p, 'T', self.TC-self.dTSC, self.medium)
        mdot = self.junctions['inlet_A'].get_massflow()

        self.junctions['outlet_A'].set_values(p=self.p, h=hout, mdot=mdot)
        self.junctions['inlet_A'].set_values(p=self.p)

    def set_air_parameters(self, T_air: float = None, mdot: float = None):
        if T_air:
            self.T_air_in = T_air
        if mdot:
            self.mdot_air = mdot

    def get_function_residual(self):
        x = np.zeros(7)
        x[0] = self.TC
        x[1] = self.TAo_desuperheat
        x[2] = self.TAo_condenser
        x[3] = self.TAo_subcool
        x[4] = self.areafraction_desuperheat
        x[5] = self.areafraction_condenser
        x[6] = self.areafraction_subcool

        res = self.model(x)
        Qdot = self.junctions['inlet_A'].get_massflow() * (self.junctions['outlet_A'].get_enthalpy() - self.junctions['inlet_A'].get_enthalpy())
        res[0:6] = res[0:6]/Qdot
        return res

class Evaporator(Component):
    def __init__(self, id: str, system: object, k: iter, area: float, superheat: float, boundary_switch: bool, limit_temp: bool):
        super().__init__(id, system)
        if len(k) == 2:
            self.k = k
        else:
            raise ValueError('k has to be of length 2. len(k) = {}'.format(len(k)))

        self.area = area
        self.superheat = superheat
        self.boundary_switch = boundary_switch
        self.limit_temp = limit_temp
        self.T0 = None
        self.TSL2 = None
        self.TSLmid = None
        self.xE1 = None
        self.xE2 = None

        self.TSL_in = None

        self.junctions['inlet_B'] = None
        self.junctions['outlet_B'] = None

    def initialize(self):
        self.p = self.junctions['inlet_A'].get_pressure()
        self.TSL1 = self.junctions['inlet_B'].get_temperature()
        self.T0 = CPPSI('T', 'P', self.p, 'Q', 0, self.junctions['inlet_A'].medium)
        self.TSL2 = self.TSL1 - (self.TSL1 - self.T0) * 0.8
        self.TSLmid = (self.TSL1 + self.TSL2)/2
        self.xE1 = 0.8
        self.xE2 = 1 - self.xE1

    def model(self, x):
        TSLi = self.junctions['inlet_B'].get_temperature()
        mSL = self.junctions['inlet_B'].get_massflow()
        hRi = self.junctions['inlet_A'].get_enthalpy()
        mR = self.junctions['inlet_A'].get_massflow()
        ref = self.junctions['inlet_A'].medium
        SL = self.junctions['inlet_B'].medium

        # print('---EVAP---')
        # print(self.junctions['inlet_A'].get_temperature())

        # Boundaries for fsolve calculation to not cause the logaritmic mean temperature to generate NaN values (neg. logarithm)
        # The refrigerants oulet temperature must not be higher than the coolants inlet temperature:

        if self.boundary_switch:
            if x[0] + self.superheat > TSLi:
                x[0] = TSLi - self.superheat

            # The evaporation temperature must not be higher than the coolants outlet temperature
            if x[0] > x[1]:
                x[0] = x[1] - 1e-6
            if x[1] > x[2]:
                x[1] = x[2] + 1e-3


        # calculate material parameters
        if (x[0] < 150.) and self.limit_temp:
            cpR = CPPSI('C', 'T', 150., 'Q', 1, ref) * (x[0]/150.)  # generate a linear extrapolation for the iteration
            hRGas = CPPSI("H", "T", 150., "Q", 1, ref) * (x[0]/150.)
        else:
            cpR = CPPSI('C', 'T', x[0], 'Q', 1, ref)  # heat capacity of fully evaporated refrigerant
            hRGas = CPPSI("H", "T", x[0], "Q", 1, ref)  # enthalpy of fully evaporated refrigerant

        cpSL = CPPSI('C', 'T', (TSLi + x[1]) / 2, 'P', 1e5, SL)  # heat capacity of secondary liquid

        # Calculate the mean logarithmic temperature value for all two sections of the condenser
        LMTD = np.zeros(2)
        LMTD[0] = LMTD_calc(x[2], x[1], x[0], x[0])
        LMTD[1] = LMTD_calc(TSLi, x[2], x[0], self.superheat + x[0])

        # Formulation of the equation system as according to fsolve documentation ( 0 = ... ).
        # The equation set  and model definition is documented in the model description.
        f = np.zeros(5)

        # energy balance evaporating zone between refrigerant and sec. liquid
        f[0] = mR * (hRGas - hRi) - mSL * cpSL * (x[2] - x[1])

        # energy balance evaporating zone between refrigerant and LMTD model
        f[1] = mR * (hRGas - hRi) - self.k[0] * x[3] * self.area * LMTD[0]

        # energy balance superheating zone between refrigerant and sec. liquid
        # f[ 2 ] = mR * (hRSuperheated - hRGas) - mSL * cpSL * (TSLi - x[ 2 ])
        f[2] = mR * cpR * self.superheat - mSL * cpSL * (TSLi - x[2])

        # energy balance superheating zone between refrigerant and LMTD model
        # f[3] = mR * (hRSuperheated-hRGas) - k[1] * x[4]/100 * Atot * LMTD[1]
        f[3] = mR * cpR * self.superheat - self.k[1] * x[4] * self.area * LMTD[1]

        # area fraction balance (1 = x_evaporating + x_superheating)
        f[4] = 1 - x[3] - x[4]

        return f

    def calc(self):
        x = np.zeros(5)
        x[0] = self.T0
        x[1] = self.TSL2
        x[2] = self.TSLmid
        x[3] = self.xE1
        x[4] = self.xE2

        x = fsolve(self.model, x0=x, xtol=self.system.fun_tol)

        self.T0 = x[0]
        self.TSL2 = x[1]
        self.TSLmid = x[2]
        self.xE1 = x[3]
        self.xE2 = x[4]

        self.p = CPPSI('P', 'T', self.T0, 'Q', 1, self.junctions['inlet_A'].medium)
        Tout = self.T0 + self.superheat
        hout = CPPSI('H', 'T', Tout, 'P', self.p, self.junctions['inlet_A'].medium)
        self.junctions['outlet_A'].set_values(p=self.p, h=hout, mdot=self.junctions['inlet_A'].get_massflow())
        hSL2 = CPPSI('H', 'T', self.TSL2, 'P', 1e5, self.junctions['inlet_B'].medium)
        # self.junctions['inlet_A'].set_values(p=self.p)
        self.junctions['outlet_B'].set_values(h=hSL2, mdot=self.junctions['inlet_B'].get_massflow())

    def get_function_residual(self):
        x = np.zeros(5)
        x[0] = self.T0
        x[1] = self.TSL2
        x[2] = self.TSLmid
        x[3] = self.xE1
        x[4] = self.xE2

        # normalize the energy balance residuals
        res = self.model(x)
        Qdot = self.junctions['inlet_A'].get_massflow() * (self.junctions['outlet_A'].get_enthalpy() - self.junctions['inlet_A'].get_enthalpy())
        res[0:4] = res[0:4]/Qdot
        return res

class IHX(Component):
    def __init__(self, id: str, system: object, UA: float):
        super().__init__(id=id, system=system)
        self.UA = UA
        self.TA_in = None
        self.TA_out = None
        self.TB_in = None
        self.TB_out = None
        self.mdot = None
        self.ref = None

        self.junctions['inlet_B'] = None
        self.junctions['outlet_B'] = None

    def initialize(self):
        self.TA_in = self.junctions['inlet_A'].get_temperature()
        self.TA_out = self.TA_in - 1.

        self.TB_in = self.junctions['inlet_B'].get_temperature()
        self.TB_out = self.TB_in + 1.

    def model(self, x):
        self.TA_in = self.junctions['inlet_A'].get_temperature()
        self.mdot = self.junctions['inlet_A'].get_massflow()
        self.TB_in = self.junctions['inlet_B'].get_temperature()
        self.medium = self.junctions['inlet_A'].medium
        # TA_out = x[0]  |  TB_out = x[1]

        pA = self.junctions['inlet_A'].get_pressure()
        pB = self.junctions['inlet_B'].get_pressure()

        cp_A = CPPSI("CPMASS", "T", self.TA_in, "P", pA, self.medium)
        cp_B = CPPSI("CPMASS", "T", self.TB_in, "P", pB, self.medium)

        if self.TA_in > self.TB_in:
            Thi = self.TA_in
            Tho = x[0]
            Tci = self.TB_in
            Tco = x[1]
        else:
            Thi = self.TB_in
            Tho = x[1]
            Tci = self.TA_in
            Tco = x[0]

        LMTD = LMTD_calc(Thi, Tho, Tci, Tco)
        Qdot = self.UA * LMTD

        f = np.zeros(2)
        if self.TA_in > self.TB_in:
            f[0] = self.mdot * cp_A * (self.TA_in - x[0]) - Qdot
            f[1] = self.mdot * cp_B * (self.TB_in - x[1]) + Qdot
        else:
            f[0] = self.mdot * cp_A * (self.TA_in - x[0]) + Qdot
            f[1] = self.mdot * cp_B * (self.TB_in - x[1]) - Qdot

        return f

    def calc(self):
        x = np.zeros(2)
        x[0] = self.TA_out
        x[1] = self.TB_out

        x = fsolve(self.model, x0=x)

        self.TA_out = x[0]
        self.TB_out = x[1]

        pA_out = self.junctions['inlet_A'].get_pressure()
        pB_out = self.junctions['inlet_B'].get_pressure()
        hA_out = CPPSI('H', 'T', self.TA_out, 'P', pA_out, self.medium)
        hB_out = CPPSI('H', 'T', self.TB_out, 'P', pB_out, self.medium)
        mdot = self.junctions['inlet_A'].get_massflow()

        self.junctions['outlet_A'].set_values(p=pA_out, h=hA_out, mdot=mdot)
        self.junctions['outlet_B'].set_values(p=pB_out, h=hB_out, mdot=mdot)
        # print('---IHX---')
        print(self.TA_out, self.TB_out)
        # print(self.junctions['outlet_A'].get_temperature(), self.junctions['outlet_B'].get_temperature())

    def get_function_residual(self):
        x = np.zeros(2)
        x[0] = self.TA_out
        x[1] = self.TB_out
        return self.model(x)


class Source(Component):
    def __init__(self, id: str, system: object, mdot=None, p=None, h=None):
        super().__init__(id=id, system=system)
        if mdot:
            self.mdot = mdot
        if p:
            self.p = p
        if h:
            self.h = h

        self.junctions = {'outlet_A': None}

    def initialize(self):
        pass

    def calc(self):
        self.junctions['outlet_A'].set_values(mdot=self.mdot, p=self.p, h=self.h)


class Sink(Component):
    def __init__(self, id: str, system: object, mdot=None, p=None, h=None):
        super().__init__(id=id, system=system)
        if mdot:
            self.mdot = mdot
        if p:
            self.p = p
        if h:
            self.h = h

        self.junctions = {'inlet_A': None}

    def initialize(self):
        pass

    def calc(self):
        pass


class HeatExchanger(Component):
    def __init__(self, id: str, system: object, UA: float):
        super().__init__(id=id, system=system)
        self.UA = UA

        self.mdotA = None
        self.TA_i = None
        self.TA_o = None
        self.TB_i = None
        self.TB_o = None
        self.mdotB = None
        self.mediumA = None
        self.mediumB = None

        self.cpA = None
        self.cpB = None

        self.junctions['inlet_B'] = None
        self.junctions['outlet_B'] = None

    def initialize(self):
        self.TA_i = self.junctions['inlet_A'].get_temperature()
        self.mdotA = self.junctions['inlet_A'].get_massflow()
        self.TB_i = self.junctions['inlet_B'].get_temperature()
        self.mdotB = self.junctions['inlet_B'].get_massflow()

        self.TA_o = self.TB_i
        self.TB_o = self.TA_i

        self.mediumA = self.junctions['inlet_A'].medium
        self.mediumB = self.junctions['inlet_B'].medium

    def model(self, x):
        # x = [TA_o, TB_o]
        # check for the hot side temperature
        if self.TA_i > self.TB_i:
            Thi = self.TA_i
            Tho = x[0]
            Tci = self.TB_i
            Tco = x[1]
        else:
            Thi = self.TB_i
            Tho = x[1]
            Tci = self.TA_i
            Tco = x[0]

        LMTD = LMTD_calc(Thi, Tho, Tci, Tco)
        Qdot = self.UA * LMTD

        f = np.zeros(2)
        if self.TA_i > self.TB_i:
            f[0] = self.mdotA * self.cpA * (self.TA_i - x[0]) - Qdot
            f[1] = self.mdotB * self.cpB * (self.TB_i - x[1]) + Qdot

        else:
            f[0] = self.mdotA * self.cpA * (self.TA_i - x[0]) + Qdot
            f[1] = self.mdotB * self.cpB * (self.TB_i - x[1]) - Qdot

        return f

    def calc(self):
        self.TA_i = self.junctions['inlet_A'].get_temperature()
        self.mdotA = self.junctions['inlet_A'].get_massflow()
        self.TB_i = self.junctions['inlet_B'].get_temperature()
        self.mdotB = self.junctions['inlet_B'].get_massflow()
        self.pA = self.junctions['inlet_A'].get_pressure()
        self.pB = self.junctions['inlet_B'].get_pressure()

        self.cpA = CPPSI('CPMASS', 'T', self.TA_i, 'P', self.pA, self.mediumA)
        self.cpB = CPPSI('CPMASS', 'T', self.TB_i, 'P', self.pB, self.mediumB)

        x = np.zeros(2)
        x[0] = self.TA_o
        x[1] = self.TB_o

        x = fsolve(self.model, x0=x)

        self.TA_o = x[0]
        self.TB_o = x[1]

        hA_o = CPPSI('H', 'T', self.TA_o, 'P', self.pA, self.mediumA)
        hB_o = CPPSI('H', 'T', self.TB_o, 'P', self.pB, self.mediumB)

        self.junctions['outlet_A'].set_values(h=hA_o)
        self.junctions['outlet_B'].set_values(h=hB_o)

        return

    def get_function_residual(self):
        x = np.zeros(2)
        x[0] = self.TA_o
        x[1] = self.TB_o
        return self.model(x)
