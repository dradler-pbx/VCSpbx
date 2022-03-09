import VCSpbx as vcs
from CoolProp.CoolProp import PropsSI as CPPSI

# parameters
SL = 'INCOMP::MEG[0.5]'
air = 'AIR'
mdot_SL = 0.25  # kg/s
p_SL = 1e5  # Pa
h_in_sl = CPPSI('H', 'T', -5. + 273.15, 'P', p_SL, SL)
mdot_air = 0.5  # kg/s
p_air = 1e5  # Pa
h_in_air = CPPSI('H', 'T', 4. + 273.15, 'P', p_air, air)

# components
system = vcs.System(id='system', tolerance=1e-4)
hex = vcs.HeatExchanger(id='hex', system=system, UA=100.)
srcSL = vcs.Source(id='srcSL', system=system, mdot=mdot_SL, p=p_SL, h=h_in_sl)
snkSL = vcs.Sink(id='snkSL', system=system)
srcAir = vcs.Source(id='srcAir', system=system, mdot=mdot_air, p=p_air, h=h_in_air)
snkAir = vcs.Sink(id='snkAir', system=system)

# junctions
srcSL_hex = vcs.Junction(id='srcSL_hex', system=system, medium=SL, upstream_component=srcSL, upstream_id='outlet_A', downstream_component=hex, downstream_id='inlet_A', mdot_init=mdot_SL, p_init=p_SL, h_init=h_in_sl)
hex_snkSL = vcs.Junction(id='hex_snkSL', system=system, medium=SL, upstream_component=hex, upstream_id='outlet_A', downstream_component=snkSL, downstream_id='inlet_A', mdot_init=mdot_SL, p_init=p_SL, h_init=h_in_sl)
srcAIR_hex = vcs.Junction(id='srcAIR_hex', system=system, medium=air, upstream_component=srcAir, upstream_id='outlet_A', downstream_component=hex, downstream_id='inlet_B', mdot_init=mdot_air, p_init=p_air, h_init=h_in_air)
hex_snkAir = vcs.Junction(id='hex_snkAir', system=system, medium=air, upstream_component=hex, upstream_id='outlet_B', downstream_component=snkAir, downstream_id='inlet_A', mdot_init=mdot_air, p_init=p_air, h_init=h_in_air)

system.run(full_output=True)
