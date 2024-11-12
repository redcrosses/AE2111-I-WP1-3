import numpy as np
from unit_conversion import *
def class_II_weight(W_dg, Nz, S_w, t_c_root, A, lam, L_t, cos_Lambda, S_csw, K_uh, F_w_Bh, S_ht, cos_Lambda_ht, A_ht, S_e, H_h, S_vt, K_2, cos_Lambda_vt, A_v, K_door, K_Lg, S_f, K_ws, L_D, K_mp, W_l, Nl, L_m, N_mw, N_mss, V_stall, K_np, N_nw, K_ng, N_Lt, W_ec, S_n, N_en,W_en, V_i,L_n, N_w, K_y,L):
    W_dg = convert_units(W_dg, 'kg', False)
    S_w = convert_units(S_w, 'm^2', False)
    L_t = convert_units(L_t, 'm', False)
    S_csw = convert_units(S_csw, 'm^2', False)
    S_ht = convert_units(S_ht, 'm^2', False)
    S_vt = convert_units(S_vt, 'm^2', False)
    S_f =convert_units(S_f, 'm^2', False)
    W_l = convert_units(W_l, 'kg', False)
    L_m = convert_units(L_m, 'm', False)
    V_stall = convert_units(V_stall, 'meters_per_sec', False)
    W_ec = convert_units(W_ec, 'kg', False)
    S_n = convert_units(S_n, 'm^2', False)
    W_en = convert_units(W_en, 'kg', False)
    V_i = convert_units(V_i, 'm^3', False)
    L_n = convert_units(L_n, 'm', False)
    L = convert_units(L, 'm', False)
    K_y = convert_units(K_y, 'm', False)




    # Wing weight calculation
    W_wing = 0.0051 * (W_dg * Nz)**0.557 * S_w**0.649 * A**0.5 * (t_c_root)**-0.4 * (1 + lam)**0.1 * (cos_Lambda)**-1 * S_csw**0.1

    # Horizontal tail weight calculation
    W_horizontal_tail = 0.0379 * K_uh * (1 + F_w_Bh)**-0.25 * W_dg**0.639 * Nz**0.10 * S_ht**0.75 * L_t**-1 * K_y**0.704 * (cos_Lambda_ht)**-1 * A_ht**0.166 * (1 + S_e / S_ht)**0.1

    # Vertical tail weight calculation
    W_vertical_tail = 0.0026 * (1 + H_h)**0.225 * W_dg**0.556 * Nz**0.536 * L_t**-0.5 * S_vt**0.5 * K_2**0.875 * (cos_Lambda_vt)**-1 * A_v**0.35 * (t_c_root)**-0.5

    # Fuselage weight calculation
    W_fuselage = 0.3280 * K_door * K_Lg * (W_dg * Nz)**0.5 * L**0.25 * S_f**0.302 * (1 + K_ws)**0.04 * (L_D)**0.10

    # Main landing gear weight calculation
    W_main_gear = 0.0106 * K_mp * W_l**0.888 * Nl**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**0.5 * V_stall**0.1

    # Nose landing gear weight calculation
    W_nose_gear = 0.032 * K_np * W_l**0.646 * Nz**0.2 * Nl**0.5 * L_n**0.5 * N_nw**0.45

    # Nacelle group weight calculation
    W_nacelle = 0.6724 * K_ng * N_Lt**0.10 * N_w **0.294 * Nz**0.611 * W_ec**0.984 * S_n**0.224

    # Starter weight calculation
    W_starter = 49.19 * ((N_en * W_en) / 1000)**0.541

    # Fuel system weight calculation
    W_fuel_system = 2.405 * V_i**0.606 * (1 + 1)**-1 * (1 + 1)**0.5

    W_uaw = 1400
    W_avionics = 1.73 * W_uaw**0.983

    W_tot = W_wing + W_horizontal_tail + W_vertical_tail + W_fuselage + W_main_gear + W_nose_gear + W_nacelle + W_starter + W_fuel_system + W_avionics

    W_tot = convert_units(W_tot, 'pounds', to_metric=True)
    W_wing = convert_units(W_wing, 'pounds', to_metric=True)
    W_fuselage = convert_units(W_fuselage, 'pounds', to_metric=True)

    return W_tot, W_wing, W_fuselage
# W_dg, Nz, S_w, t_c_root, A, lam, L_t, cos_Lambda, S_csw, K_uh, F_w_Bh, S_ht, cos_Lambda_ht, A_ht, S_e, H_h, S_vt, K_2, cos_Lambda_vt, A_v, K_door, K_Lg, S_f, K_ws, L_D, K_mp, W_l, Nl, L_m, N_mw, N_mss, V_stall, K_np, N_nw, K_ng, N_Lt, W_ec, S_n, N_en,W_en, V_i,L_n, N_w, K_y,L = 110389.61038961037, 3.75, 134.81210140752708, 0.12, 8.059564745668995, 0.30090020507176196 ,16.3, 0.9999626055152955, 91.0692369597575, 1, 0, 33.14589928240281, 0.7880107536067219, 4, 0, 0, 19.139135064849572, 16.3, 0.766044443118978, 1.5, 1.06, 1, 287.0108363467119, 0,7732328503748613, 15.70558015983471, 1, 96038.96103896103, 4.5, 6, 12, 2, 69.44, 1, 2, 1.017, 4.5, 5900, 0, 2, 5900, 1, 6, 0, 4.89, 21.276940375

# W_tot, W_wing, W_fuselage = class_II_weight(W_dg, Nz, S_w, t_c_root, A, lam, L_t, cos_Lambda, S_csw, K_uh, F_w_Bh, S_ht, cos_Lambda_ht, A_ht, S_e, H_h, S_vt, K_2, cos_Lambda_vt, A_v, K_door, K_Lg, S_f, K_ws, L_D, K_mp, W_l, Nl, L_m, N_mw, N_mss, V_stall, K_np, N_nw, K_ng, N_Lt, W_ec, S_n, N_en,W_en, V_i,L_n, N_w, K_y,L)

# Variables with unique names and their definitions

# W_dg = W_dg           # Design gross weight
# Nz = N_z            # Ultimate load factor (1.5 × limit load factor)
# S_w = S_optimal          # Wing area
# t_c_root = t_cratio      # Thickness-to-chord ratio at root
# lam = taper_ratio          # Wing taper ratio
# cos_Lambda = np.cos(np.radians(sweep_quarter))    # Cosine of wing sweep angle at 25% MAC
# S_csw = S_wf         # Control surface area (wing-mounted)
# K_y =0.3*L_t
# K_uh = 1           # Empirical factor for horizontal tail weight
# F_w_Bh = 0        # Wing-fuselage interference factor and horizontal tail span factor
# S_ht = htail_area          # Horizontal tail area
# cos_Lambda_ht = np.cos(radians(htail_sweep)) # Cosine of horizontal tail sweep angle
# A_ht =  htail_AR         # Aspect ratio of horizontal tail
# S_e = 0           # Elevator area

# H_h = 0            # Empirical factor for vertical tail weight
# S_vt = vtail_area          # Vertical tail area
# K_2 = L_t           # Vertical tail volume coefficient
# cos_Lambda_vt = np.cos(radians(vtail_sweep)) # Cosine of vertical tail sweep angle
# A_v = vtail_AR           # Aspect ratio of vertical tail

# K_door = 1.06        # Door factor for fuselage weight
# K_Lg = 1          # Gear location factor for fuselage weight
# S_f =  S_wfuselage          # Fuselage wetted area
# K_ws = 0.75*(1+2*taper_ratio)/(1+taper_ratio)(span*tan(sweep_quarter)/(l_cyl))          # Wing sweep factor
# L_D = liftoverdrag           # Lift-to-drag ratio

# K_mp = 1          # Main gear weight factor
# W_l = Wl           # Landing design gross weight
# L_m = L_m           # Main gear wheelbase
# N_mw = 12          # Number of main wheels
# N_mss = 2         # Number of main gear shock struts
# V_stall = stall_speed       # Stall speed

# K_np = 1          # Nose gear weight factor
# N_ln = L_n          # Length of nose landing gear
# N_nw = 2          # Number of nose wheels

# K_ng = 1.017          # Nacelle weight factor
# N_Lt = Nl          # Ultimate landing load factor (1.5 × N_gear)
# W_ec = 5900          # Weight of engine and contents (per nacelle)
# S_n = 0           # Nacelle wetted area

# Nl = Nl
# W_en = 5900          # Engine weight
# V_i = volume_f           # Integral tanks volume


