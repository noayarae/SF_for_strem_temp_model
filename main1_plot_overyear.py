### This code compute the Shade-Factor for streams

### Input
# - "in01_topog_angle.csv":        CSV file with topographic angles for 0, 45, 90, 135, 270, 225, 270, 315 axis
# - "in02_stream_data_LM_RM.csv":  CSV file with stream data: (0)sb,(1)long_rch,(2)lat_rch, (3)treeh_r, (4)treeh_l, (5)strm_az, (6)strm_w, (7)strm_slen
# - "in03_SR_data_2019.csv":       CSV file with daily Solar Radiation: (0)Date, (1)julian_d, (2)SR
# - Set time step "t_step" at L.150. It is suggested to set 0.01
# - Set riparian "gap" at L.126      It is suggested to set 10

### Output
# - "sf001_top.csv"
# - "sf002_Rtree.csv"
# - "sf003_Ltree.csv"
# - "sf004_all.csv"

### SF over the year is plotted calling another function

# import required modules

import csv, math,time#, sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
#from random import randrange

import get_sun_az_angle as sa_ang  # call own functions/files

start_time = datetime.now()


### Read atmospheric angle data
l_atm_ang = []
with open('in01_topog_angle.csv', 'r') as file:
    # (0)0 deg, (1)45, (2)90, (3)135, (4)180, (5)225, (6)270, (7)315
    reader1 = csv.reader(file)
    next(reader1)
    for row in reader1:
        at_ang1 = map(float, list(row))
        at_ang2 = at_ang1[1:]
        l_atm_ang.append(at_ang2)
        #print(at_ang2)

    #print len(l_atm_ang)
    #print (l_atm_ang)
### --------- End: Read atm data -----------------

#print ssdd

### Read stream data
strm_db = []

### with open('in02_stream_data_for_sf.csv', 'r') as file:
# (0) with open('in02_stream_data_L0_R0.csv', 'r') as file:
# (1) with open('in02_stream_data_LN_RN.csv', 'r') as file:
# (2) with open('in02_stream_data_LM_RM.csv', 'r') as file:
# (3) with open('in02_stream_data_LR_effic.csv', 'r') as file:
# (4) with open('in02_stream_data_LR_effic2.csv', 'r') as file:
sf_case = ['in02_stream_data_L0_R0.csv','in02_stream_data_LN_RN.csv','in02_stream_data_LM_RM.csv',
           'in02_stream_data_LR_effic.csv','in02_stream_data_LR_effic2.csv']
sf_id = 3  #  <-------------- SET

#with open('in02_stream_data_LR_effic.csv', 'r') as file:
#with open('in02_stream_data_LR_effic.csv', 'r') as file:
with open(sf_case[sf_id], 'r') as file:
    # (1)Long,(2)Lat, (3)R_h-tree, (4)L_h-tree, (5)Str_Az, (6)Str_W, (7)Str_Len
    reader2 = csv.reader(file)
    next(reader2)
    #OutputList = [map(float, list(i)) for i in zip(*reader)]
    for row in reader2:
        strm_d1 = map(float, list(row))
        strm_d2 = strm_d1[1:]
        strm_db.append(strm_d2)
        #print(strm_d2)
        
    ### (0)Long,(1)Lat, (2)R_h-tree, (3)L_h-tree, (4)Str_Az, (5)Str_W, (6)Str_Len
    #print len(strm_db)
    #print (strm_db)
### --------- End: Read stream data -----------------

#print asdas

### Read SR data
sr_db = []
with open('in03_SR_data_2019.csv', 'r') as file: # Solar Radiation Data
    # (1)date, (2)SR
    reader3 = csv.reader(file)
    next(reader3)
    #OutputList = [map(float, list(i)) for i in zip(*reader)]
    for row in reader3:
        sr1 = row[1:]
        sr2 = map(float, list(sr1))
        sr_db.append(sr2[1])
        #print(sr2)

    print "L57 SR data length:",len(sr_db)
    #print (sr_db)
### --------- End: Read SR data -----------------

#print kk


### Time information - Days # 1/1/2019: 43466   12/31/2019: 43830
#j_day = 40350 # Day: 6/21/2010
#t_past = 0.2/24. # time
day001 = 43466#43466#43637 #43466 # Jan-01-2019
#day172 = 43637; print "6/21/2019" # Day: 6/21/2019
day365 = 43830 #43830 # Dec-31-2019
#days = day365-day001
lt_days = np.arange(day001, day365+1, 1) # list_time_days
lt_n_day = np.arange(1, (day365-day001+2), 1) # list_number_of_time_days
print "Range of days: [",day001,"-",day365,"]", " Number_of_days =", len(lt_days)
#print lt_days
#print lt_n_day

sf_topf = []
sf_Lf = []
sf_Rf = []
sf_allf = []


print "len-stream DB:",len(strm_db)
### (0)Long,(1)Lat, (2)R_h-tree, (3)L_h-tree, (4)Str_Az, (5)Str_W, (6)Str_Len

n_sbs = 60#60#len(strm_db)
### For loop through SBs
print range(n_sbs)
#print dasdasd
for sb in [20]:#range(n_sbs):#[28]:#  #  <-----  It may SET
    print sb  # It starts from 0 to 59

    ### Site information
    longg = strm_db[sb][0] # -123.206# -122.65 # -105 # Longitude
    latt = strm_db[sb][1] # 45.734# 45.5169444 # 40 # Latitude
    t_zone = -8 # -6  # Time zone

    ### Stream vegetation data
    gap = 5. #10.                 # <--- SET
    az_str = strm_db[sb][4] # 90. # <--- SET
    hr_tree = strm_db[sb][2]# 10. # <--- SET
    hl_tree = strm_db[sb][3]# 10. # <--- SET

    sf_yr_top = []
    sf_yr_R_tree = []
    sf_yr_L_tree = []
    sf_yr_RL_tree = []
    sf_yr_all = []

    ### For loop through Days
    for dy in range(len(lt_days)):
        #print "dy:",dy
        if (lt_n_day[dy]%20 ==0 or lt_n_day[dy]==365): print "SB",(sb+1)," Day",lt_n_day[dy], "  Julian_d:", lt_days[dy]
        #print "Day",lt_n_day[dy], "   Julian_day:", lt_days[dy]
        pot_dy_sr = 0 # Potential day Solar Rad
        sr_dy_blk_top = 0 # SR day blocked by Topography
        sr_dy_blk_L_tree = 0 # SR day blocked by L-tree
        sr_dy_blk_R_tree = 0 # SR day blocked by R-tree
        sr_dy_blk_RL_tree = 0 # SR day blocked by RL-tree
        sr_dy_blk_all = 0 # SR day blocked by All
        #print adsa

        ### Time information - Hours and fractions
        t_step = 0.01#0.01#0.1 #1 # 0.1 means fraction of 1 hours  # <-------- SET
        #print "Time step:", t_step
        lt_past = np.arange(0, 24, t_step)
        #print len(lt_past)
        #print "Number of day-fractions: ",lt_past

        list_sa = []
        list_az = []
        #print "Solar Angle"
        #print "Solar Azimuth"
            
        for tp in lt_past: # Loop through each fraction of time within 1-day
            #print "L.110 tp:",tp
            #listd = list_input = [latt, longg, t_zone, j_day, tp/24.]
            listd = list_input = [latt, longg, t_zone, lt_days[dy], tp/24.]
            sa,az,sunset = sa_ang.calc_sun_az_angle (listd)  # ---->  Call function
            list_sa.append(sa)
            list_az.append(az)
            #print "L.109 sa & sunset:",sa, sunset
            #print az

            ### Calc hourly SolRad
            tss = 24.*sunset # day = lt_n_day[dy]
            term1 = 24.*lt_n_day[dy]* math.sin(2*math.pi*tss/24)/math.pi
            term2 = 2.*lt_n_day[dy]*(tss-12)* math.cos(2*math.pi*tss/24)
            bb = sr_db[dy]/(term1 - term2)
            aa = -bb*math.cos(2*math.pi*tss/24)

            # Solar radiation at tp instant
            sr_hr = lt_n_day[dy] * aa + bb * lt_n_day[dy] * math.cos(2*math.pi*tp/24)
            #print "test"
            #print tp, bb, aa, sr_hr
            

            if sa>0: # If solar angle > 0 do calculations
                
                #if az>360: print("Alert, Azimuth angle > 360") # Alert if Azimuth is > 360
                #qdrt = int((az+22.5-360)/45) if (az+22.5>=360) else int((az+22.5)/45) # Get Quadrant
                #print "Sun_ang: ", sa, "Az: ", az, "  Az-sector", qdrt, "  top_ang", l_atm_ang[0][qdrt]

                ### Get top angle for the azimuth by interpolation
                #az = 360.
                if az == 360: az = 359.99 # Makes a bit lower than 360 to avoid crash
                idx1 = int(az/45)
                idx2 = idx1-7 if (idx1+1>7) else idx1+1
                #print "uppper, lower bound:", l_atm_ang[0][idx1], l_atm_ang[0][idx2] # upper and lower bound of top_angle
                ratio = (l_atm_ang[sb][idx2]-l_atm_ang[sb][idx1])/45.
                top_ang = l_atm_ang[sb][idx1] + (az-idx1*45.)*ratio
                #print az, ratio, top_ang
                #print "Sun_ang: ", sa, "Az: ", az, "  top_ang", top_ang

                #print wwee
                
                ### Potential hourly SolRad
                if sr_hr<-1:
                    print "Alert SR negative but Sol-Angle posite"
                    print "L.151 stop:",stop
                elif sr_hr>-1 and sr_hr < 0:
                    sr_hr = 0
                else:
                    pot_hr_sr = strm_db[sb][5] *1.0* sr_hr #pot_hr_sr = strm_db[0][5] *1.0* sr_hr
                    #print "pot_sr:", strm_db[0][5], sr_hr, pot_hr_sr
                    pot_dy_sr += pot_hr_sr
                ### End Potential hourly SolRad

                #if sa < l_atm_ang[0][qdrt]: # Solar-angle < topograph-angle
                if sa < top_ang: # Solar-angle < topograph-angle
                    #print "Do calc",dy, strm_db[0][5], sr_hr
                    #shade_top = strm_db[0][5]*strm_db[0][6]*sr_db[dy]
                    SR_blk_top = strm_db[sb][5] * 1.0 * sr_hr #SR_blk_top = strm_db[0][5] * 1.0 * sr_hr # Width*Len*SR
                    #print "SR_blocked by topog:",SR_blk_top
                    sr_dy_blk_top += SR_blk_top
                    #print sds
                else:            # Solar-angle > topog-angle
                    if math.sin((az-az_str)*math.pi/180)>0: # R-side
                        sd1 = hr_tree/math.tan(sa*math.pi/180)
                    else:                                   # L-side
                        sd1 = hl_tree/math.tan(sa*math.pi/180)
        
                    #sd1 = hr_tree/math.tan(sa*math.pi/180)
                    sd2 = sd1*math.sin((az-az_str)*math.pi/180)
                    #print "sd1, sd2:", sd1, sd2

                    #if abs(sd2)>(strm_db[0][5]+gap): # Shade >> stream-width. Stream is fully shaded
                    if abs(sd2)>(strm_db[sb][5]+gap): # Shade >> stream-width. Stream is fully shaded
                        if sd2 < 0:     # Shade is caused by Left-Trees
                            #print "case --- 1"
                            #SR_blk_L = strm_db[0][5] * 1.0 * sr_hr # Width*Len*SR
                            SR_blk_L = strm_db[sb][5] * 1.0 * sr_hr # Width*Len*SR
                            SR_blk_R = 0
                        else:           # Shade is caused by Right-Trees
                            #print "case --- 2"
                            SR_blk_L = 0
                            #SR_blk_R = strm_db[0][5] * 1.0 * sr_hr # Width*Len*SR
                            SR_blk_R = strm_db[sb][5] * 1.0 * sr_hr # Width*Len*SR
                    elif abs(sd2)< gap:              # Shade < Gap.           Stream is no-shaded
                        #print "case --- 3" # No shade over the stream because the shade is lower than gap distance
                        SR_blk_L = 0
                        SR_blk_R = 0
                    else:                            # Gap < Shade < Width.  Stream is partially shaded
                        if sd2 < 0:     # Shade is caused by Left-Trees
                            #print "case --- 4"
                            SR_blk_L = (abs(sd2)-gap) * 1.0 * sr_hr #
                            SR_blk_R = 0
                        else:           # Shade is caused by Right-Trees
                            #print "case --- 5"
                            SR_blk_L = 0
                            SR_blk_R = (abs(sd2)-gap) * 1.0 * sr_hr #
                            
                    #print "SR_blocked by tree R,L:",SR_blk_L, SR_blk_R
                    sr_dy_blk_L_tree += SR_blk_L
                    sr_dy_blk_R_tree += SR_blk_R
                    sr_dy_blk_RL_tree += SR_blk_L + SR_blk_R
                    
        sr_dy_blk_all = sr_dy_blk_top + sr_dy_blk_RL_tree
                        
        #print "pot_day_Sol-Rad:", pot_dy_sr
        #print "SR blocked by topog:", sr_dy_blk_top # SR day blocked by topog
        #print "SR blocked by R-tree:", sr_dy_blk_R_tree # SR day blocked by R-tree
        #print "SR blocked by L-tree:", sr_dy_blk_L_tree # SR day blocked by L-tree
        #print "SR blocked by RL-tree:", sr_dy_blk_RL_tree # SR day blocked by RL
        #print "SR blocked by All:", sr_dy_blk_all # SR day blocked by All

        sf_top = sr_dy_blk_top/pot_dy_sr
        sf_R_tree = sr_dy_blk_R_tree/pot_dy_sr
        sf_L_tree = sr_dy_blk_L_tree/pot_dy_sr
        sf_RL_tree = sr_dy_blk_RL_tree/pot_dy_sr
        sf_all = sr_dy_blk_all/pot_dy_sr
        #print "SF by topog:", sf_top # SF day by topog
        #print "SF by R-tree:", sf_R_tree # SF day by R-tree
        #print "SF by L-tree:", sf_L_tree # SF day by L-tree
        #print "SF by RL-tree:", sf_RL_tree # SF day by RL
        #print "SF by All:", sf_all # SF day BY All
                    
        #print sds
            
        #print "Number of time fractions: ",len(list_sa)# Numero de fracciones de calculo
        #print "list of solar angles (Neg & Pos)"
        #print list_sa # list of solar angles including negatives
        #print np.sign(list_sa) # funcion signo
        #print list_az

        ### Get sub-set of sun-ang and az-ang for only positives
        range_idx = []
        for idx in range(0, len(list_sa) - 1):
            # checking for successive opposite index
            if list_sa[idx] > 0 and list_sa[idx+1] < 0 or list_sa[idx] < 0 and list_sa[idx+1] > 0:
                range_idx.append(idx)
        #print "Index of range of positive sun_angles: ", range_idx

        #subset_sa = list_sa[range_idx[0]+1:range_idx[1]+1]
        #subset_az = list_az[range_idx[0]+1:range_idx[1]+1]
        #print "Len of Pos Sol_Ang:", len(subset_sa)
        #print "List of Positive solar angles"
        #print subset_sa
        #print "List of azimuths for positive Sol_Ang"
        #print subset_az
        #print

        sf_yr_top.append(sf_top)
        sf_yr_R_tree.append(sf_R_tree)
        sf_yr_L_tree.append(sf_L_tree)
        sf_yr_RL_tree.append(sf_RL_tree)
        sf_yr_all.append(sf_all)
        

        '''
        ### ---- Plot Solar radiation over the day ---- ###
        import matplotlib.pyplot as plt
        xpoints = np.array(lt_past)
        ypoints = np.array(list_sa)

        plt.plot(xpoints, ypoints)
        plt.show()
        #'''

    from pprint import pprint
    print "Shade Factor by Topographic"
    #print sf_yr_top# pprint (sf_yr_top)
    print "Shade Factor by R-Tree"
    #print sf_yr_R_tree
    print "Shade Factor by L-Tree"
    #print sf_yr_L_tree
    print "Shade Factor by RL-Tree"
    #print sf_yr_RL_tree
    print "Shade Factor by ALL"
    #print sf_yr_all

    ### set Shade-Factor final files
    sf_topf.append(sf_yr_top)
    sf_Rf.append(sf_yr_R_tree)
    sf_Lf.append(sf_yr_L_tree)
    sf_allf.append(sf_yr_all)

    ### Comment these two lines to turn-off plot
    import plot_areas_overyear as py #  ------------- FUNCTION  PLOT
    plot_ar = py.plot_areas (sb, sf_yr_top, sf_yr_R_tree, sf_yr_L_tree, sf_yr_all, strm_db[sb][4],sf_id)
    '''
    ### ---- Plot Shade factor over the year ---- ###
    ### L.156 t_step = 0.1
    ### Data to plot
    x_pts = np.array(list(range(1,366)))
    top_pts = np.array(sf_yr_top)
    Rtree_pts = np.array(sf_yr_R_tree)
    Ltree_pts = np.array(sf_yr_L_tree)
    all_pts = np.array(sf_yr_all)

    ### Basic setings of the plot
    fig = plt.figure(figsize=(8, 4), facecolor='white') # Figure size, background color
    ### Setings of the 4 lines: Top, Rtree, Ltree, All
    p1, = plt.plot(x_pts, top_pts, color='b', label='Top',dashes=[4,2], linestyle = ':', linewidth = '2')
    p2, = plt.plot(x_pts, Rtree_pts, color='g', label='R-tree',dashes=[10,3], linestyle = '--', linewidth = '2')
    p3, = plt.plot(x_pts, Ltree_pts, color='r', label='L-tree',dashes=[12,2,3,2,3,2], linestyle = '-.', linewidth = '2')
    p4, = plt.plot(x_pts, all_pts, color='k', label='All', linestyle = '-', linewidth = '2')
    ### Text for legend
    txt1 = "Shade factor by topographic"
    txt2 = "Shade-factor by Right bank riparian vegetation"
    txt3 = "Shade-factor by Left bank riparian vegetation"
    txt4 = "Shade-factor by all obstacles"
    leg = plt.legend([p1, p2, p3, p4], [txt1, txt2,txt3,txt4], loc='upper left',
                     scatterpoints=1, prop={'size': 9}, bbox_to_anchor=(0.01,0.999),handlelength=4.5)

    plt.xticks(np.arange(0, 370, 20))
    plt.xlim([0, 370])
    plt.ylim([0, 1.42])
    plt.grid()

    title = 'Shade-factor by different obstacles for sub-basin '+str(sb+1)
    fig.suptitle(title, fontsize=12)
    plt.xlabel('Day in the year', fontsize=10)
    plt.ylabel('Shade Factor (SF)', fontsize=10)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=10)
    leg.get_lines()[0].set_linewidth(2.5)
    leg.get_lines()[1].set_linewidth(2.5)
    leg.get_lines()[2].set_linewidth(2.5)
    leg.get_lines()[3].set_linewidth(2.5)    
    
    plt.show()
    #'''
print


sf_topft = map(list, zip(*sf_topf))  # Transpose of "sf_topf"
sf_Rft = map(list, zip(*sf_Rf))  # Transpose of "sf_R"
sf_Lft = map(list, zip(*sf_Lf))  # Transpose of "sf_L"
sf_allft = map(list, zip(*sf_allf))  # Transpose of "sf_allf"
print "N_days & N_SBs: ",len(sf_topft), len(sf_topft[0])

head = []
for i in range(len(sf_topft[0])):
    head.append("sb"+str(i+1))
#print head

### Write Files for top-SF, R-SF, L-SF
with open("sf001_top.csv", "wb") as f1:
    writer = csv.writer(f1)
    #writer.writerows([["sb1","sb2"]])
    writer.writerows([head])
    writer.writerows(sf_topft)

with open("sf002_Rtree.csv", "wb") as f2:
    writer = csv.writer(f2)
    #writer.writerows([["sb1","sb2"]])
    writer.writerows([head])
    writer.writerows(sf_Rft)

with open("sf003_Ltree.csv", "wb") as f3:
    writer = csv.writer(f3)
    #writer.writerows([["sb1","sb2"]])
    writer.writerows([head])
    writer.writerows(sf_Lft)

with open("sf004_all.csv", "wb") as f4:
    writer = csv.writer(f4)
    #writer.writerows([["sb1","sb2"]])
    writer.writerows([head])
    writer.writerows(sf_allft)

#'''

print "done..."
end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))












