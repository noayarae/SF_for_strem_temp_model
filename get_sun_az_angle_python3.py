### This python function computes the Solar azimute angle
# It can be run by itself.
# This is called by the "main1_SF_calculator_overyear_python3.py"

import math
#import csv, math, sys
#import numpy as np
#from random import randrange
#import matplotlib.pyplot as plt

def input_data():
    ### Site information
    latt = 45.5169444 #45.6145477 #40 # Latitude
    longg = -122.65 #-123.0715942 #-105 # Longitude
    t_zone = -8 #-6  # Time zone

    ### Time information
    j_day = 40350 # Day: 6/21/2010
    #j_day = 43466; print "1/1/2019" # Day: 1/1/2019
    j_day = 43637; print ("6/21/2019") # Day: 6/21/2019
    #j_day = 43830; print "12/31/2019" # Day: 12/31/2019

    t_past = 8/24. # time
    print ("t_past: ", t_past)
    list_input = [latt, longg, t_zone, j_day, t_past]
    return list_input

def r2d(ang_rad):
    ang_deg = ang_rad * 180./math.pi
    return ang_deg
def d2r(ang_deg):
    ang_rad = ang_deg*math.pi/180.
    return ang_rad

def calc_sun_az_angle (list_input):
    latt = list_input[0] # Latitude
    longg = list_input[1] # Longitude
    t_zone = list_input[2]  # Time zone
    day = list_input[3] # Day
    t_past = list_input[4] # time
    
    # Julian day, Julian century
    julian_d = day + 2415018.5 + t_past - t_zone/24.
    julian_c = (julian_d-2451545)/36525
    long_s = (280.46646 + julian_c*(36000.76983 + julian_c*0.0003032))%360
    anom_s = 357.52911 + julian_c*(35999.05029 - 0.0001537*julian_c)
    #print "julian_d: ",julian_d
    #print "julian_c: ",julian_c
    #print "long_s: ",long_s
    #print "anom_s: ",anom_s

    # Get eccent Earth orbit
    excent = 0.016708634 - julian_c*(0.000042037 + 0.0000001267*julian_c)

    sun_eq1 = math.sin(d2r(anom_s))*(1.914602-julian_c*(0.004817+0.000014*julian_c))
    sun_eq2 = math.sin(d2r(2*anom_s))*(0.019993-0.000101*julian_c)
    sun_eq3 = math.sin(d2r(3*anom_s))*0.000289
    sun_eq = sun_eq1 + sun_eq2 + sun_eq3
    sun_tlong = long_s + sun_eq

    sun_tanom = anom_s + sun_eq
    sun_radv = (1.000001018*(1-excent*excent))/(1+excent*math.cos(d2r(sun_tanom)))
    sun_along = sun_tlong-0.00569-0.00478*math.sin(d2r(125.04-1934.136*julian_c))

    #print "excent: ",excent
    #print "sun_eq: ",sun_eq
    #print "sun_tlong: ",sun_tlong
    #print "sun_tanom: ", sun_tanom
    #print "sun_radv: ", sun_radv
    #print "sun_along: ", sun_along

    # Get Mean obliq ecliptic (deg) and others
    obl_e = 23+(26+((21.448-julian_c*(46.815+julian_c*(0.00059-julian_c*0.001813))))/60)/60
    obl_cor = obl_e + 0.00256*math.cos(d2r(125.04-1934.136*julian_c))
    sun_asc1 = math.cos(d2r(sun_along))
    sun_asc2 = math.cos(d2r(obl_cor))*math.sin(d2r(sun_along))
    sun_asc = r2d(math.atan2(sun_asc2,sun_asc1))

    sun_dec = r2d(math.asin(math.sin(d2r(obl_cor))*math.sin(d2r(sun_along))))
    var_y = math.tan(d2r(obl_cor/2))*math.tan(d2r(obl_cor/2))

    #print "obl_e: ", obl_e
    #print "obl_cor: ", obl_cor
    #print "sun_asc: ", sun_asc
    #print "sun_dec: ", sun_dec
    #print "var_y: ", var_y

    # Get et_min and others
    et_min1 = var_y*math.sin(2*d2r(long_s)) - 2*excent*math.sin(d2r(anom_s))
    et_min2 = 4*excent*var_y*math.sin(d2r(anom_s))*math.cos(2*d2r(long_s))
    et_min3 = -0.5*var_y*var_y*math.sin(4*d2r(long_s)) - 1.25*excent*excent*math.sin(2*d2r(anom_s))
    et_min = 4*r2d(et_min1 + et_min2 + et_min3)

    sun_rise1 = math.cos(d2r(90.833))/(math.cos(d2r(latt))*math.cos(d2r(sun_dec)))
    sun_rise2 = math.tan(d2r(latt))*math.tan(d2r(sun_dec))
    sun_rise = r2d(math.acos(sun_rise1-sun_rise2))

    sol_noon = (720 - 4*longg - et_min + t_zone*60)/1440
    sunrise_t = sol_noon - sun_rise*4/1440.
    sunset_t = sol_noon + sun_rise*4/1440.
    sun_dur = 8*sun_rise

    #print "et_min: ", et_min
    #print "sun_rise: ", sun_rise
    #print "sol_noon: ", sol_noon
    #print "sunrise_t: ", sunrise_t
    #print "sunset_t: ", sunset_t
    #print "sun_dur: ", sun_dur

    # Get true_st and others
    true_st = (t_past*1440 + et_min + 4*longg - 60*t_zone) % 1440

    #IF(AB3/4<0,AB3/4+180,AB3/4-180)
    if true_st/4.<0:
      #print("op_1")
      hour_ang = true_st/4. + 180
    else:
      #print("op_2")
      hour_ang = true_st/4. - 180


    sol_zang1 = math.sin(d2r(latt))*math.sin(d2r(sun_dec))
    sol_zang2 = math.cos(d2r(latt))*math.cos(d2r(sun_dec))*math.cos(d2r(hour_ang))
    sol_zang = r2d(math.acos(sol_zang1 + sol_zang2))

    sol_elev = 90 - sol_zang

    #print "true_st: ", true_st
    #print "hour_ang: ", hour_ang
    #print "hour_ang: ", sol_zang
    #print "sol_elev: ", sol_elev

    ### Get at_refr and others
    #sol_elev = 85.1 # To test
    if sol_elev > 85:
      #print("op_1")
      at_ref = 0
    else:
      #print("op_2")
      if sol_elev > 5:
          #print("op_2a")
          at_ref1 = 58.1/math.tan(d2r(sol_elev))-0.07/math.pow(math.tan(d2r(sol_elev)),3)
          at_ref2 = 0.000086/math.pow(math.tan(d2r(sol_elev)),5)
          at_ref = at_ref1 + at_ref2
      else:
          #print("op_2b")
          if sol_elev > -0.575:
              #print("op_2b-1")
              at_ref = 1735 + sol_elev*(-518.2 + sol_elev*(103.4 + sol_elev*(-12.79 + sol_elev*0.711)))
          else:
              #print("op_2b-2")
              at_ref = -20.772/math.tan(d2r(sol_elev))

    at_ref = at_ref/3600.
    #print "sol_elev: ", sol_elev
    #print "at_ref: ", at_ref

    ### Get sol_elev_cor and others
    sol_elev_cor = sol_elev + at_ref

    saa1 = ((math.sin(d2r(latt))*math.cos(d2r(sol_zang))) - math.sin(d2r(sun_dec)))
    saa2 = (math.cos(d2r(latt))*math.sin(d2r(sol_zang)))
    if hour_ang > 0:
        #print("ch_1")
        sol_az_ang = (r2d(math.acos(saa1/saa2))+180) % 360
    else:
        #print("ch_2")
        sol_az_ang = (540-r2d(math.acos(saa1/saa2))) % 360

    #print "sol_elev_cor: ", sol_elev_cor
    #print "sol_az_ang: ", sol_az_ang

    return sol_elev_cor, sol_az_ang, sunset_t

if __name__ == "__main__":
    # Input data
    list_input = input_data()

    sun_angle, az_angle, sunset = calc_sun_az_angle(list_input)
    print ("Solar-angle ; Azimut-angle:")
    print (sun_angle, az_angle)

    print ("donee...")
















