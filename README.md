### Code to compute the Shade-Factor (SF) 

This code computes the riparian SF, which was added to the Ficklin stream temperature model (I called it the Improved Ficklin stream temperature model) and later to include into the SWAT model v2012 (rev681). This approach has been used in analyzes of Effects of shade-factor from riparian vegetation in the Dairy-McKay Basin, Oregon, USA on stream water temperature, which were shown in the publication "An improved model of shade-affected stream temperature in Soil & Water Assessment Tool"
https://doi.org/10.5194/hess-27-739-2023

Input data to run it are:

- "in01_topog_angle.csv":        CSV file with topographic angles for azimuth: 0, 45, 90, 135, 270, 225, 270, 315. Zero for the North
- "in01_topog_angle.csv":        CSV file with topographic angles for the azimute: 0, 45, 90, 135, 270, 225, 270, 315
- "in02_stream_data_LM_RM.csv":  CSV file with stream data: (0)number_of_sb,(1)longitude_rch,(2)latitude_rch, (3)tree_height_right_bank, (4)tree_hight_left_bank, (5)strm_azimute, (6)strm_width, (7)strm_slength
- "in03_SR_data_2019.csv":       CSV file with daily Solar Radiation: (0)Date, (1)julian_day, (2)Solar_Radiation

Some of the following settings are suggested
- Set time step "t_step" at around L.167. It is suggested to set 0.01
- Set riparian "gap" at around L.142      Here, I set 5 meters
- Set 'list_sbs' at aroud L.131:          with the list of sub-basins. If you want to run this for all the sun-basins, set a list with their number. For instance,
  * if you watershed has 60 sub-basins, set:    list_sbs = list(range(0, n_sbs)) 
  * if you want to run only for sub-basins 21 and 26, set: list_sbs = [20, 25] 
