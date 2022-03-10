import VCSpbx as vcs
from CoolProp.CoolProp import PropsSI as CPPSI
from fmipp.export.createFMU import createFMU

# parameter setting
cpr_speed = 2400.0
T_amb = 30. + 273.15
ref = 'R290'
sl = 'INCOMP::MEG[0.5]'
h_in_sl = CPPSI('H', 'P', 1e5, 'T', 273.15 - 1.85, sl)
mdot_SL = 0.25
p_SL = 1e5
superheat = 4.0

mdot_init = 4e-3
pc_init = 14e5
Tc_init = CPPSI('T', 'P', pc_init, 'Q', 0, ref)
h2_init = CPPSI('H', 'P', pc_init, 'T', 60+273.15, ref)
h3_init = CPPSI('H', 'P', pc_init, 'Q', 0, ref)
h4_init = CPPSI('H', 'P', pc_init, 'T', Tc_init-2., ref)
p0_init = 2.5e5
T0_init = CPPSI('T', 'P', p0_init, 'Q', 1, ref)
h5_init = CPPSI('H', 'P', p0_init, 'T', T0_init+superheat, ref)
h1_init = CPPSI('H', 'P', p0_init, 'T', T0_init+superheat+2., ref)
h_out_SL_init = CPPSI('H', 'P', 1e5, 'T', 273.15 - 10.0, sl)

system = vcs.System(id='system', tolerance=1e-4)
cpr = vcs.Compressor_efficiency(id='cpr', system=system, etaS=0.65, etaV=0.9, stroke=33e-6, speed=cpr_speed)
cond = vcs.Condenser(id='cond', system=system, k=[450., 450., 450.], area=1., subcooling=0.1, T_air_in=T_amb, mdot_air_in=0.56)
ihx = vcs.IHX(id='ihx', system=system, UA=2.3)
evap = vcs.Evaporator(id='evap', system=system, k=[420., 420.], area=1., superheat=superheat, boundary_switch=True, limit_temp=True)
srcSL = vcs.Source(id='srcSL', system=system, mdot=mdot_SL, p=p_SL, h=h_in_sl)
snkSL = vcs.Sink(id='snkSL', system=system)

cpr_cond = vcs.Junction(id='cpr_cond', system=system, medium=ref, upstream_component=cpr, upstream_id='outlet_A', downstream_component=cond, downstream_id='inlet_A', mdot_init=mdot_init, p_init=pc_init, h_init=h2_init)
cond_ihx = vcs.Junction(id='cond_ihx', system=system, medium=ref, upstream_component=cond, upstream_id='outlet_A', downstream_component=ihx, downstream_id='inlet_A', mdot_init= mdot_init, p_init=pc_init, h_init=h3_init)
ihx_evap = vcs.Junction(id='ihx_evap', system=system, medium=ref, upstream_component=ihx, upstream_id='outlet_A', downstream_component=evap, downstream_id='inlet_A', mdot_init= mdot_init, p_init=p0_init, h_init=h4_init)
evap_ihx = vcs.Junction(id='evap_ihx', system=system, medium=ref, upstream_component=evap, upstream_id='outlet_A', downstream_component=ihx, downstream_id='inlet_B', mdot_init=mdot_init, p_init=p0_init, h_init=h5_init)
ihx_cpr = vcs.Junction(id='ihx_cpr', system=system, medium=ref, upstream_component=ihx, upstream_id='outlet_B', downstream_component=cpr, downstream_id='inlet_A', mdot_init=mdot_init, p_init=p0_init, h_init=h1_init)
srcSL_evap = vcs.Junction(id='srcSL_evap', system=system, medium=sl, upstream_component=srcSL, upstream_id='outlet_A', downstream_component=evap, downstream_id='inlet_B', mdot_init=mdot_SL, p_init=p_SL, h_init=h_in_sl)
evap_snkSL = vcs.Junction(id='evap_snkSL', system=system, medium=sl, upstream_component=evap, upstream_id='outlet_B', downstream_component=snkSL, downstream_id='inlet_A', mdot_init=mdot_SL, p_init=p_SL, h_init=h_out_SL_init)

# system.run(full_output=True)

fmu_name = 'test_vcs_FMU'

fmu_dict = {'Inputs':
                {'n_cpr': cpr.set_speed,
                 'T_amb': cond.set_air_temp,
                 'mdot_cond': cond.set_air_mdot,
                 'h_SL': srcSL.set_enthalpy,
                 'mdot_SL': srcSL.set_mdot},
            'Outputs':
                {'T_aircond_out': cond.get_outlet_temp,
                 'T_SL_out': snkSL.get_temperature,
                 'P_el': cpr.get_power}}
