from django.shortcuts import render
import math
from django.http import HttpResponse

#Functions

def surface_data(u):
    

    if u == 1:
        return [9.54, 13.95, 2.68, 0.254, 1250, 0.840]
    elif u == 2:
        return [6.48, 16.00, 1.86, 0.152, 1804, 0.845]
    elif u == 3:
        return [5.21, 19.82, 1.54, 0.102, 2231, 0.841]
    elif u == 4:
        return [5.11, 20.06, 1.49, 0.102, 2290, 0.843]
    else:
        print("Wrong Info")
        return []

def alpha(B1, B2, beta, a=0.1524):
    cal = B1 + B2 + 2 * a
    ans = B1 / cal
    return ans * beta

def sigma(alpha, hydraulic_dia):
    return (hydraulic_dia * alpha) / (4 * 1000)

def reynold_no(mass, actual_area, hydraulic_dia, d_viscosity):
    return (hydraulic_dia * mass * 10) / (actual_area * d_viscosity)

def prandtl_no(d_viscosity, specific_heat, conductivity):
    return (d_viscosity * specific_heat) / (conductivity * 10000)

def heat_transfer_coef(colburn_fact, mass, actual_area, specific_heat, prandtl_no):
    cal = math.pow(prandtl_no, 2.0 / 3.0)
    return (colburn_fact * specific_heat * mass) / (actual_area * cal)

def fin_eff(heat_transfer_coeff, metal_conduct, thickness):
    cal = (2.0 * heat_transfer_coeff * 1000) / (metal_conduct * thickness)
    cal = math.pow(cal, 0.5)
    return cal

def net_fin_efficiency(fin_eff, fin_length):
    u = 2 * fin_eff * fin_length
    o = math.exp(u)
    return (o - 1) / ((o + 1) * (u / 2))

def overall_fluid_eff(net_fin_eff, surface_ratio):
    cal = (1 - net_fin_eff)
    cal *= surface_ratio
    return (1 - cal)

def overall_heat_tc(overall_eff_c, surface_area_c, heat_transfer_coef_c, overall_eff_h, surface_area_h, heat_transfer_coef_h):
    cal = (1 / (overall_eff_c * surface_area_c * heat_transfer_coef_c))
    cal += (1 / (overall_eff_h * surface_area_h * heat_transfer_coef_h))
    return (1 / (surface_area_c * cal))

def n_t_u(overall_eff, surface_area, min_specific_heat, mass):
    return (overall_eff * surface_area) / (min_specific_heat * mass)

def loss_1(Kc, sigma):
    return 1 + Kc - sigma * sigma

def loss_2(volume_out, volume_in):
    return 2 * ((volume_out / volume_in) - 1)

def loss_3(friction_fact, surface_area, actual_area, volume_out, volume_in):
    return friction_fact * surface_area * ((volume_out / volume_in) + 1) / (2 * actual_area)

def loss_4(ke, sigma, volume_out, volume_in):
    return (1 - sigma * sigma - ke) * volume_out / volume_in





def home_view(request):


    context = {}

    # Take form values
    if request.method == 'POST':
        user_data = request.POST
        # print(user_data)
        
        Temp_c_in = float(user_data['Temp_c_in'])
        mass_c = float(user_data['mass_c'])
        specific_heat_c = float(user_data['specific_heat_c'])*1000
        d_viscos_c = float(user_data['d_viscos_c'])
        conduct_c = float(user_data['conduct_c'])
        Temp_h_in = float(user_data['Temp_h_in'])
        Temp_h_out = float(user_data['Temp_h_out'])
        mass_h = float(user_data['mass_h'])
        specific_heat_h = float(user_data['specific_heat_h'])*1000
        d_viscos_h = float(user_data['d_viscos_h'])
        conduct_h = float(user_data['conduct_h'])
        height = float(user_data['height'])
        depth = float(user_data['depth'])
        width = float(user_data['width'])
        op = int(user_data['op'])  # Assuming op is an integer
        metal_conduct = float(user_data['metal_conduct'])
        press_in_c = float(user_data['press_in_c'])
        press_drop_all_c = float(user_data['press_drop_all_c'])
        press_in_h = float(user_data['press_in_h'])
        press_drop_all_h = float(user_data['press_drop_all_h'])




        #Calculations

        surface_data_c_value = int(request.POST.get('surface_data_c', None))
        surface_data_h_value = int(request.POST.get('surface_data_h', None))

        # print(surface_data_c_value)
        # print(surface_data_h_value)


        surface_data_c = []
        surface_data_h = []

        if surface_data_c_value is not None:
            surface_data_c = surface_data(surface_data_c_value)
            # Now you can use the surface_data_c variable as needed
        else:
            print("Cold surface data value is not provided or is None.")

        if surface_data_h_value is not None:
            surface_data_h = surface_data(surface_data_h_value)
        else:
            print("Hot surface data value is not provided or is None.")
 
        Heat = mass_h * ((Temp_h_in - Temp_h_out) / 1000.0) * specific_heat_h
        Temp_c_out = Temp_c_in + (Heat * 1000 / (mass_c * specific_heat_c))

        volume = height * depth * width

        prandtl_no_c = prandtl_no(d_viscos_c, specific_heat_c, conduct_c)
        prandtl_no_h = prandtl_no(d_viscos_h, specific_heat_h, conduct_h)

        if op  == 1:
            frontal_area_c = width * height
            frontal_area_h = depth * height
        else:
            frontal_area_h = width * height
            frontal_area_c = depth * height

        alpha_c = alpha(surface_data_c[0], surface_data_h[0], surface_data_c[4])
        alpha_h = alpha(surface_data_h[0], surface_data_c[0], surface_data_h[4])
        
        total_surface_area_c = alpha_c * volume
        total_surface_area_h = alpha_h * volume
        
        sigma_c = sigma(alpha_c, surface_data_c[2])
        sigma_h = sigma(alpha_h, surface_data_h[2])
        
        actual_area_c = frontal_area_c * sigma_c
        actual_area_h = frontal_area_h * sigma_h


        reynold_no_c = reynold_no(mass_c, actual_area_c, surface_data_c[2], d_viscos_c)
        reynold_no_h = reynold_no(mass_h, actual_area_h, surface_data_h[2], d_viscos_h)
        
        j_factor_c = 0.47 / math.pow(reynold_no_c, 0.5)
        j_factor_h = 0.47 / math.pow(reynold_no_h, 0.5)


        kec, kc, keh, kh, fc, fh = 0, 0, 0, 0, 0, 0
        
        if reynold_no_c > 2000:
            kc = 0.6 - 0.45 * sigma_c
            kec = 1 - 1.15 * sigma_c
            fc = 0.26 / math.pow(reynold_no_c, 0.25)
        else:
            kc = 1.3 - 0.35 * sigma_c
            kec = 1 - 1.9 * sigma_c
            fc = 0.43 / math.pow(reynold_no_c, 0.3)
        
        if reynold_no_h > 2000:
            kh = 0.6 - 0.45 * sigma_h
            keh = 1 - 1.15 * sigma_h
            fh = 0.26 / math.pow(reynold_no_h, 0.25)
        else:
            kh = 1.3 - 0.35 * sigma_h
            keh = 1 - 1.9 * sigma_h
            fh = 0.43 / math.pow(reynold_no_h, 0.3)


        h_t_c = heat_transfer_coef(j_factor_c, mass_c, actual_area_c, specific_heat_c, prandtl_no_c)
        h_t_h = heat_transfer_coef(j_factor_h, mass_h, actual_area_h, specific_heat_h, prandtl_no_h)
        
        
        fin_eff_c = fin_eff(h_t_c, metal_conduct, surface_data_c[3])
        fin_eff_h = fin_eff(h_t_h, metal_conduct, surface_data_h[3])
        
        leng_fin = 3.175 / 1000.0
        
        net_fin_eff_c = net_fin_efficiency(fin_eff_c, leng_fin)
        overall_fluid_eff_c = overall_fluid_eff(net_fin_eff_c, surface_data_c[5])
        net_fin_eff_h = net_fin_efficiency(fin_eff_h, leng_fin)
        overall_fluid_eff_h = overall_fluid_eff(net_fin_eff_h, surface_data_h[5])
        
        overall_heat_coef, NTU, effectiveness = 0, 0, 0
        
        if mass_c * specific_heat_c > mass_h * specific_heat_h:
            overall_heat_coef = overall_heat_tc(overall_fluid_eff_h, total_surface_area_h, h_t_h, overall_fluid_eff_c, total_surface_area_c, h_t_c)
            NTU = n_t_u(overall_heat_coef, total_surface_area_h, specific_heat_h, mass_h)
            effectiveness = (Temp_h_in - Temp_h_out) / (Temp_h_in - Temp_c_in)
        else:
            overall_heat_coef = overall_heat_tc(overall_fluid_eff_c, total_surface_area_c, h_t_c, overall_fluid_eff_h, total_surface_area_h, h_t_h)
            NTU = n_t_u(overall_heat_coef, total_surface_area_c, specific_heat_c, mass_c)
            effectiveness = (Temp_c_out - Temp_c_in) / (Temp_h_in - Temp_c_in)

        press_out_c = press_in_c - press_drop_all_c
        
        input_vol_c = (287 * (Temp_c_in + 273)) / (press_in_c * 1000)
        output_vol_c = (287 * (Temp_c_out + 273)) / (press_out_c * 1000)
        
        cal = loss_1(kc, sigma_c) + loss_2(output_vol_c, input_vol_c) + loss_3(fc, total_surface_area_c, actual_area_c, output_vol_c, input_vol_c) - loss_4(kec, sigma_c, output_vol_c, input_vol_c)
        press_loss_c = (mass_c * mass_c) / (actual_area_c * actual_area_c)
        press_loss_c *= input_vol_c / (2 * 9.80665)
        press_loss_c *= cal
        press_loss_c /= 1000


        press_out_h = press_in_h - press_drop_all_h


        input_vol_h = (287 * (Temp_h_in + 273)) / (press_in_h * 1000)
        output_vol_h = (287 * (Temp_h_out + 273)) / (press_out_h * 1000)
        
        cal2 = loss_1(kh, sigma_h) + loss_2(output_vol_h, input_vol_h) + loss_3(fh, total_surface_area_h, actual_area_h, output_vol_h, input_vol_h) - loss_4(keh, sigma_h, output_vol_h, input_vol_h)
        press_loss_h = (mass_h * mass_h) / (actual_area_h * actual_area_h)
        
        press_loss_h *= input_vol_h / (2 * 9.80665)
        press_loss_h /= 1000
        press_loss_h *= cal2


        # print("Total Expected Heat Exchange between fluid (Q) (in KW): ", Heat)
        # print("Exit Temperature of the Cold Fluid is (in *C): ", Temp_c_out)
        # print("Pressure loss for the Cold fluid (in Kpa): ", press_loss_c)
        # print("Pressure loss for the Hot fluid (in Kpa): ", press_loss_h)
        # print("NTU of system: ", NTU)
        # print("Overall Heat transfer coefficient of the system is (in W/m²k): ", overall_heat_coef)





        context = {
            'Heat': round(Heat, 3),
            'Temp_c_out': round(Temp_c_out, 3),
            'press_loss_c': round(press_loss_c, 3),
            'press_loss_h': round(press_loss_h, 3),
            'NTU': round(NTU, 3),
            'overall_heat_coef': round(overall_heat_coef, 3),
            'effectiveness': effectiveness
        }










    return render(request, 'index.html', context)