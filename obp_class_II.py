from numpy import *
from unit_conversion import *

class class_II_weight:
  def __init__(self, W_dg, Nz, S_w, t_c_root, A, lam, L_t, cos_Lambda, S_csw, K_uh, F_w_Bh, S_ht, cos_Lambda_ht, A_ht, S_e, H_h, S_vt, K_2, cos_Lambda_vt, A_v, K_door, K_Lg, S_f, K_ws, L_D, K_mp, W_l, Nl, L_m, N_mw, N_mss, V_stall, K_np, N_nw, K_ng, N_Lt, W_ec, S_n, N_en,W_en, V_i,L_n, N_w, K_y,L, powerplant_mass) -> None:
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
    V_i = V_i * 264.172
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
    W_main_gear = 0.0106 * K_mp * W_l**0.888 * Nl**0.25 * L_m**0.4 * N_mw**0.321 * N_mss**-0.5 * V_stall**0.1

    # Nose landing gear weight calculation
    W_nose_gear = 0.032 * K_np * W_l**0.646 * Nz**0.2 * Nl**0.5 * L_n**0.5 * N_nw**0.45

    # Nacelle group weight calculation
    W_nacelle = 0.6724 * K_ng * N_Lt**0.10 * N_w **0.294 * Nz**0.611 * W_ec**0.984 * S_n**0.224

    # Starter weight calculation
    W_starter = 49.19 * ((N_en * W_en) / 1000)**0.541

    # Fuel system weight calculation
    W_fuel_system = 2.405 * V_i**0.606 * (1 + 1)**-1 * (1 + 1)**0.5 * 3

    W_uaw = 1400
    W_avionics = 1.73 * W_uaw**0.983

    #W_tot = W_wing + W_horizontal_tail + W_vertical_tail + W_fuselage + W_main_gear + W_nose_gear + W_nacelle + W_starter + W_fuel_system + W_avionics




    self.wing = convert_units(W_wing,'pounds', True)
    self.htail = convert_units(W_horizontal_tail,'pounds', True)
    self.vtail = convert_units(W_vertical_tail,'pounds', True)
    self.empennage = self.htail+self.vtail
    self.fus = convert_units(W_fuselage,'pounds', True)
    self.main_lg = convert_units(W_main_gear,'pounds', True)
    self.nose_lg = convert_units(W_nose_gear,'pounds', True)
    self.gear = self.nose_lg + self.main_lg
    self.nacelle = convert_units(W_nacelle,'pounds', True)
    self.starter = convert_units(W_starter,'pounds', True)
    self.fuel_sys = convert_units(W_fuel_system,'pounds', True)
    self.powerplant = powerplant_mass
    self.systems = self.powerplant + self.fuel_sys

    self.total = 1.33 * (self.wing + self.htail +self.vtail + self.fus + self.main_lg + self.nose_lg + self.fuel_sys + self.powerplant) #added a proportion to account for removable parts; proportion taken from ADSEE II lecture 5 slides #+ self.nacelle + self.starter
    pass

  def printall(self):
    print("wing:",self.wing, "rest:",self.htail,self.vtail,self.fus,self.main_lg,self.nose_lg,self.fuel_sys, self.powerplant)
    
# values = [110389.61038961037, 3.75, 134.81210140752708, 0.12, 8.059564745668995, 0.30090020507176196, 16.3, 0.9999626055152955, 91.0692369597575, 1, 0, 33.14589928240281, 0.7880107536067219, 4, 0, 0, 19.139135064849572, 16.3, 0.766044443118978, 1.5, 1.06, 1, 287.0108363467119, 0.7732328503748613, 15.70558015983471, 1, 96038.96103896103, 4.5, 6, 12, 2, 69.44, 1, 2, 1.017, 4.5, 5900, 0, 2, 5900, 58.36850649350649, 6, 0, 4.89, 21.276940375]

# (W_dg, Nz, S_w, t_c_root, A, lam, L_t, cos_Lambda, S_csw, K_uh, F_w_Bh, S_ht, cos_Lambda_ht, A_ht, S_e, H_h, S_vt, K_2, cos_Lambda_vt, A_v, K_door, K_Lg, S_f, K_ws, L_D, K_mp, W_l, Nl, L_m, N_mw, N_mss, V_stall, K_np, N_nw, K_ng, N_Lt, W_ec, S_n, N_en, W_en, V_i, L_n, N_w, K_y, L) = values

# W_class_II = class_II_weight(W_dg, Nz, S_w, t_c_root, A, lam, L_t, cos_Lambda, S_csw, K_uh, F_w_Bh, S_ht, cos_Lambda_ht, A_ht, S_e, H_h, S_vt, K_2, cos_Lambda_vt, A_v, K_door, K_Lg, S_f, K_ws, L_D, K_mp, W_l, Nl, L_m, N_mw, N_mss, V_stall, K_np, N_nw, K_ng, N_Lt, W_ec, S_n, N_en,W_en, V_i,L_n, N_w, K_y,L)

# print("wing:", W_class_II.wing)
# print("total:", W_class_II.total)
# print(lam)
