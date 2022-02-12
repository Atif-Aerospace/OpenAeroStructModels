"""
Conversion of related geometry parameters in different ways of wing defintion
0 = given absolute values of span and chords
1 = given wing area and all the ratios
2 = Trapezoid wing, no kink
"""

import math

# 0 = given absolute values of span and chords
def SizeRelation_0(entire_span, kink_location, root_chord, kink_chord, tip_chord):
    
    # wing entire span in m
    entirespan_dependency = "(Independent)" 
    # note if this varaible is given or computed, due to multiple options to define the wing

    # wing semi-span in m
    half_span = entire_span / 2
    halfspan_dependency = "(Dependent)"

    # spanwise location of the kink from the center line in m
    kink_location_ratio = kink_location / half_span
    kinklocation_dependency = "(Independent)"
    kinklocation_ratio_dependency = "(Dependent)"

    # root chord in m
    rootchord_dependency = "(Independent)"

    # kink chord in m
    kinkchord_dependency = "(Independent)"

    # tip chord in m
    tipchord_dependency = "(Independent)"

    # Taper ratio at tip/root
    taper_ratio_t_r = tip_chord / root_chord
    taperratioTr_dependency = "(Dependent)"

    # Taper ratio at kink/root
    taper_ratio_k_r = kink_chord / root_chord
    taperratioKr_dependency = "(Dependent)"

    # wing area in m^2
    wing_area_geo = (root_chord + kink_chord) * kink_location + (kink_chord + tip_chord) * (half_span - kink_location)
    wingarea_geo_dependency = "(Dependent)"

    # wing aspect ratio
    aspect_ratio_geo = entire_span ** 2 / wing_area_geo
    aspectratio_geo_dependency = "(Dependent)"

    return (half_span, kink_location_ratio, taper_ratio_t_r, taper_ratio_k_r, wing_area_geo, aspect_ratio_geo, 
    wingarea_geo_dependency, aspectratio_geo_dependency, entirespan_dependency, halfspan_dependency, kinklocation_dependency, kinklocation_ratio_dependency,
    taperratioTr_dependency, taperratioKr_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency)

# 1 = given wing area and all the ratios
def SizeRelation_1(wing_area_geo, aspect_ratio_geo, kink_location_ratio, taper_ratio_t_r, taper_ratio_k_r):
    # wing area in m^2
    wingarea_geo_dependency = "(Independent)"

    # wing aspect ratio
    aspectratio_geo_dependency = "(Independent)"

    # wing entire span in m
    entire_span = math.sqrt(wing_area_geo * aspect_ratio_geo)
    entirespan_dependency = "(Dependent)" 

    # wing semi-span in m
    half_span = entire_span / 2
    halfspan_dependency = "(Dependent)"

    # spanwise location of the kink from the center line in m
    kink_location = kink_location_ratio * half_span
    kinklocation_ratio_dependency = "(Independent)"
    kinklocation_dependency = "(Dependent)"

    # Taper ratio at tip/root
    taperratioTr_dependency = "(Independent)"

    # Taper ratio at kink/root
    taperratioKr_dependency = "(Independent)"

    # root chord in m
    root_chord = wing_area_geo / (kink_location * (1 + taper_ratio_k_r) + (half_span - kink_location) * (taper_ratio_k_r + taper_ratio_t_r))
    rootchord_dependency = "(Dependent)"

    # kink chord in m
    kink_chord = root_chord * taper_ratio_k_r
    kinkchord_dependency = "(Dependent)"

    # tip chord in m
    tip_chord = root_chord * taper_ratio_t_r
    tipchord_dependency = "(Dependent)"

    return (entire_span, half_span, kink_location, root_chord, kink_chord, tip_chord, 
    wingarea_geo_dependency, aspectratio_geo_dependency, entirespan_dependency, halfspan_dependency, kinklocation_dependency, kinklocation_ratio_dependency,
    taperratioTr_dependency, taperratioKr_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency)

# 2 = Trapezoid wing, no kink
def SizeRelation_2(wing_area_geo, aspect_ratio_geo, taper_ratio_t_r):
    
    # wing area in m^2
    wingarea_geo_dependency = "(Independent)"

    # wing aspect ratio
    aspectratio_geo_dependency = "(Independent)"

    # wing entire span in m
    entire_span = math.sqrt(wing_area_geo * aspect_ratio_geo)
    entirespan_dependency = "(Dependent)" 

    # wing semi-span in m
    half_span = entire_span / 2
    halfspan_dependency = "(Dependent)"

    # spanwise location of the kink from the center line in m
    kink_location = half_span
    kink_location_ratio = 1
    kinklocation_dependency = "(Dependent)"
    kinklocation_ratio_dependency = "(Dependent)"

    # Taper ratio at tip/root
    taperratioTr_dependency = "(Independent)"

    # Taper ratio at kink/root
    taper_ratio_k_r = taper_ratio_t_r
    taperratioKr_dependency = "(Same as taper_ratio_t_r)"

    # root chord in m
    root_chord = wing_area_geo / (half_span * (1 + taper_ratio_t_r))
    rootchord_dependency = "(Dependent)"

    # kink chord in m
    kink_chord = root_chord * taper_ratio_k_r
    kinkchord_dependency = "(Same as taper_ratio_t_r)"

    # tip chord in m
    tip_chord = root_chord * taper_ratio_t_r
    tipchord_dependency = "(Dependent)"

    return (entire_span, half_span, kink_location, kink_location_ratio, taper_ratio_k_r, root_chord, kink_chord, tip_chord, 
    wingarea_geo_dependency, aspectratio_geo_dependency, entirespan_dependency, halfspan_dependency, kinklocation_dependency, kinklocation_ratio_dependency, 
    taperratioTr_dependency, taperratioKr_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency)
    

# Wimpress definition of the wing
def SizeRelation_Wimpress_1(wing_area_wpref, aspect_ratio_wpref, kink_location_ratio, body_side_ratio, taper_ratio_trap, root_chord_extension_ratio):
    
    # Wimpress reference area in m^2
    wingarea_wpref_dependency = "(Independent)"

    # Aspect ratio based on Wimpress reference area
    aspectratio_wpref_dependency = "(Independent)"

    # wing entire span in m
    entire_span = math.sqrt(wing_area_wpref * aspect_ratio_wpref)
    entirespan_dependency = "(Dependent)" 

    # wing semi-span in m
    half_span = entire_span / 2
    halfspan_dependency = "(Dependent)"

    # spanwise location of the kink from the center line in m
    kink_location = kink_location_ratio * half_span
    kinklocation_ratio_dependency = "(Independent)"
    kinklocation_dependency = "(Dependent)"

    # Root chord of trap area
    # see D:\Tools\OpenAeroStruct\openaerostruct\examples\mytest\pic\Wimpress Definition.jpg
    temp = (1 + taper_ratio_trap) + (root_chord_extension_ratio - 1) * kink_location_ratio * ((kink_location_ratio - body_side_ratio) / (1 - body_side_ratio))
    root_chord_trap = 2 * wing_area_wpref / (entire_span * temp)
    bodysideratio_dependency = "(Independent)"
    taperratioTrap_dependency = "(Independent)"
    rootchordTrap_dependency = "(Dependent)"


    # The actual root chord
    root_chord = root_chord_extension_ratio * root_chord_trap
    rootchord_extratio_dependency = "(Independent)"
    rootchord_dependency = "(Dependent)"

    # tip chord in m
    tip_chord = root_chord_trap * taper_ratio_trap
    tipchord_dependency = "(Dependent)"

    # kink chord in m
    kink_chord = (1 - kink_location_ratio) * (root_chord_trap - tip_chord) + tip_chord
    kinkchord_dependency = "(Dependent)"

    # Taper ratio at tip/root
    taper_ratio_t_r = tip_chord / root_chord
    taperratioTr_dependency = "(Dependent)"

    # Taper ratio at kink/root
    taper_ratio_k_r = kink_chord / root_chord
    taperratioKr_dependency = "(Dependent)"

    # geometry (gross) wing area in m^2
    wing_area_geo = (root_chord + kink_chord) * kink_location + (kink_chord + tip_chord) * (half_span - kink_location)
    wingarea_geo_dependency = "(Dependent)"

    # geometry (gross) wing aspect ratio
    aspect_ratio_geo = entire_span ** 2 / wing_area_geo
    aspectratio_geo_dependency = "(Dependent)"

    # trap wing area in m^2
    wing_area_trap = (root_chord_trap + tip_chord) * half_span
    wingarea_trap_dependency = "(Dependent)"

    # trap wing aspect ratio
    aspect_ratio_trap = entire_span ** 2 / wing_area_trap
    aspectratio_trap_dependency = "(Dependent)"


    return (wing_area_geo, aspect_ratio_geo, wing_area_trap, aspect_ratio_trap, entire_span, half_span, kink_location, 
    root_chord_trap, root_chord, kink_chord, tip_chord, taper_ratio_t_r, taper_ratio_k_r,
    wingarea_wpref_dependency, aspectratio_wpref_dependency, 
    wingarea_geo_dependency, aspectratio_geo_dependency, 
    wingarea_trap_dependency, aspectratio_trap_dependency,
    kinklocation_ratio_dependency, bodysideratio_dependency, rootchord_extratio_dependency, taperratioTrap_dependency, taperratioTr_dependency, taperratioKr_dependency,
    entirespan_dependency, halfspan_dependency, kinklocation_dependency, rootchordTrap_dependency, rootchord_dependency, kinkchord_dependency, tipchord_dependency)

# 0 = 0.25 chord line sweep
def SweepRelation_Wimpress_0(trap_quarter_sweep, kink_location, root_chord_trap, kink_chord, root_chord):
    # trap area quarter chord sweep angle in deg
    trapquartersweep_dependency = "(Independent)"

    # x location of the quarter point of the kink chord (from root chord leading edge)
    x_kink_quarter = 0.25 * root_chord_trap + math.tan(trap_quarter_sweep / 180.0 * math.pi) * kink_location

    # x location of the leading edge point of the kink chord (measured from root chord leading edge)
    x_kink_LE = x_kink_quarter - 0.25 * kink_chord

    # x location of the trailing edge point of the kink chord (measured from root chord leading edge)
    x_kink_TE = x_kink_quarter + 0.75 * kink_chord

    # inboard leading-edge sweep angle in deg
    inboard_LE_sweep = math.atan(x_kink_LE / kink_location) / math.pi * 180
    inLEsweep_dependency = "(Dependent)"
    outboard_LE_sweep = inboard_LE_sweep
    outLEsweep_dependency = "(Dependent)"


    # inboard quarter chord sweep angle in deg
    inquartersweep_dependency = "(Dependent)"
    inboard_quarter_sweep = math.atan( (x_kink_quarter - 0.25 * root_chord) / kink_location ) / math.pi * 180

    # outboard quarter chord sweep angle in deg
    outquartersweep_dependency = "(Dependent)"
    outboard_quarter_sweep = trap_quarter_sweep

    return (inboard_quarter_sweep, outboard_quarter_sweep, inboard_LE_sweep, outboard_LE_sweep, 
    trapquartersweep_dependency, inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency)




# 0 = 0.25 chord line sweep
def SweepRelation_0(inboard_quarter_sweep, outboard_quarter_sweep, half_span, kink_location, root_chord, kink_chord, tip_chord):
    # inboard quarter chord sweep angle in deg
    inquartersweep_dependency = "(Independent)"

    # outboard quarter chord sweep angle in deg
    outquartersweep_dependency = "(Independent)"

    # x location of the quarter point of the kink chord (from root chord leading edge)
    x_kink_quarter = 0.25 * root_chord + math.tan(inboard_quarter_sweep / 180.0 * math.pi) * kink_location

    # x location of the leading edge point of the kink chord (from root chord leading edge)
    x_kink_LE = x_kink_quarter - 0.25 * kink_chord

    # inboard leading-edge sweep angle in deg
    inboard_LE_sweep = math.atan(x_kink_LE / kink_location) / math.pi * 180
    inLEsweep_dependency = "(Dependent)"

    # x location of the quarter point of the tip chord (from root chord leading edge)
    x_tip_quarter = x_kink_quarter + math.tan(outboard_quarter_sweep / 180.0 * math.pi) * (half_span - kink_location)

    # x location of the leading edge point of the tip chord (from root chord leading edge)
    x_tip_LE = x_tip_quarter - 0.25 * tip_chord

    # outboard leading-edge sweep angle in deg
    if kink_location == half_span:
        outboard_LE_sweep = inboard_LE_sweep
        outLEsweep_dependency = "(Same as inboard_LE_sweep)"
    else:
        outboard_LE_sweep = math.atan((x_tip_LE - x_kink_LE) / (half_span - kink_location)) / math.pi * 180
        outLEsweep_dependency = "(Dependent)"
    
    return (inboard_LE_sweep, outboard_LE_sweep, 
    inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency)

# 1 = leading edge sweep
def SweepRelation_1(inboard_LE_sweep, outboard_LE_sweep, half_span, kink_location, root_chord, kink_chord, tip_chord):
    # inboard leading-edge sweep angle in deg
    inLEsweep_dependency = "(Independent)"

    # outboard leading-edge sweep angle in deg
    outLEsweep_dependency = "(Independent)"

    # x location of the leading edge point of the kink chord (from root chord leading edge)
    x_kink_LE = math.tan(inboard_LE_sweep / 180.0 * math.pi) * kink_location

    # x location of the quarter point of the kink chord (from root chord leading edge)
    x_kink_quarter = x_kink_LE + 0.25 * kink_chord

    # inboard leading-edge sweep angle in deg
    inboard_quarter_sweep = math.atan( (x_kink_quarter - 0.25 * root_chord) / kink_location ) / math.pi * 180
    inquartersweep_dependency = "(Dependent)"

    # x location of the leading edge point of the tip chord (from root chord leading edge)
    x_tip_LE = x_kink_LE + math.tan(outboard_LE_sweep / 180.0 * math.pi) * (half_span - kink_location)

    # x location of the quarter point of the tip chord (from root chord leading edge)
    x_tip_quarter = x_tip_LE + 0.25 * tip_chord

    # outboard leading-edge sweep angle in deg
    if kink_location == half_span:
        outboard_quarter_sweep = inboard_quarter_sweep
        outquartersweep_dependency = "(Same as inboard_quarter_sweep)"
    else:
        outboard_quarter_sweep = math.atan( (x_tip_quarter - x_kink_quarter) / (half_span - kink_location) ) / math.pi * 180
        outquartersweep_dependency = "(Dependent)"
    
    return (inboard_quarter_sweep, outboard_quarter_sweep, 
    inquartersweep_dependency, outquartersweep_dependency, inLEsweep_dependency, outLEsweep_dependency)