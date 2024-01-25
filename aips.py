import sys

import casacore.tables as pt
from casacore.tables import table, taql
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import AltAz, SkyCoord, FK5, EarthLocation, get_body


def get_elevation_and_azimuth(source, start_time_, stop_time_, coordinates_of_source):
    solar_system_objects = ["sun", "jupiter"]
    coordinates_rt32 = EarthLocation(x=3183649.31400 * u.m, y=1276902.98900*u.m, z=5359264.71000 * u.m)

    start_time = Time(start_time_, scale='utc')
    stop_time = Time(stop_time_, scale='utc')
    duration = int((stop_time.mjd - start_time.mjd) * 24 * 60 * 60)
    start_time = Time(start_time_, scale='utc')
    time_range = np.linspace(0, duration, duration) * u.second
    times = start_time + time_range

    if source in solar_system_objects:
        body = get_body(source, times, coordinates_of_lofar)
        body = body.transform_to(FK5)
        ra = body.ra
        dec = body.dec

    else:
        if source == "J2253+417":
            source_ = SkyCoord(ra="22h46m49.7323000s", dec="44d20'02.368000", frame=FK5, equinox='J2000.0')
        elif source == "EV_LAC":
            source_ = SkyCoord(ra="22h55m36.7078000s", dec="42d02'52.533000", frame=FK5, equinox='J2000.0')

    
    frame = AltAz(obstime=times, location=coordinates_rt32)
    elevation_azimuth = source_.transform_to(frame)

    return elevation_azimuth


def get_tsys(antab_file):
    tsys = dict()
    with open(antab_file, "r") as antab:
         antab_data = antab.readlines()
         for line in antab_data:
             if "source=" in line:
                 src = line.split("source=")[1].replace("\n", "").upper()
                 if src not in tsys:
                     
                     tsys[src] = {str(ch):[] for ch in range(0, 16)}
             elif line.startswith("!") or "INDEX=" in line or "TSYS IR FT" in line or "/" in line:
                continue
             else:
                 t_sys = line.split(" ")
                 for ch in range(2, 18):
                    tsys[src][str(ch-2)].append(np.float(t_sys[ch].replace("\n", "")))
                
    return tsys


def main():
    source_J2253_417_scans = [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73, 76]
    source_EV_LAC_scans = [2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23, 24, 26, 27, 29, 30, 32, 33, 35, 36, 38, 39, 41, 42, 44, 45, 47, 48, 50, 51, 53, 54, 56, 57, 59, 60, 62, 63, 65, 66, 68, 69, 71, 72, 74, 75]


    elevation_poly = [4.5642e-06, -0.0007491198, 1.0132377024]
    dfpu = [0.085300, 0.082000]

    msname = "stef24_contiuum.ms"
    freqs = pt.table(msname+"/SPECTRAL_WINDOW").getcol("CHAN_FREQ") /10**6
    ms_table = table(msname)
  
    field_ids = list(set(ms_table.getcol("FIELD_ID")))
    field_table = msname + "/FIELD"
    field = table(field_table)
    sources = field.getcol("NAME")
    coordinates = field.getcol("PHASE_DIR")

    print("sources: ", sources)

    time = ms_table.getcol("TIME")/60/60/24
    time2 = Time(np.array(ms_table.getcol('TIME'))/60/60/24, format='mjd', scale='utc')
    time_datetime = time2.to_datetime()
    time_datetime = [dt.strftime("%Y_%B_%d_%H:%M:%S") for dt in time_datetime]
    
    start_time = Time(np.min(time2), format='mjd', scale='utc').to_datetime().strftime("%Y-%m-%dT%H:%M:%S") 
    stop_time = Time(np.max(time2), format='mjd', scale='utc').to_datetime().strftime("%Y-%m-%dT%H:%M:%S") 
    
    print("start time: ", np.min(time2), start_time)
    print("stop time:", np.max(time2), stop_time)
    
    nr_if = len(freqs)
    frequencies_list_for_each_if = [np.linspace(np.min(freqs[i]), np.max(freqs[i]), len(freqs[i])) for i in range(0, nr_if)]
    print("Number of IF: ", nr_if)

    antab_file = msname.split("_")[0] + "ir.antabfs"
    print("antab file", antab_file)
    tsys = get_tsys(antab_file)
    tsys_ch_map = {0:[1,9],1:[2,10],2:[3,11],3:[4,12],4:[5,13],5:[6,14],6:[7,15],7:[8,16]}
    
    output_calibrator = []
    output_calibrator_poly1d = []
    output_target = []
    calibrator = False 
    for source in sources:
        if source == "J2253+417":
            calibrator = True            

        print(source)
        source_index = sources.index(source)
        source_table = taql('select from $ms_table where FIELD_ID == ' + str(source_index))
        data_for_source_raw = source_table.getcol("DATA")
        print(data_for_source_raw.shape)

        flag = source_table.getcol("FLAG")
        data_for_source_raw = np.ma.array(data_for_source_raw, mask=flag)
      
        coords = coordinates[source_index][0]
        elevation = get_elevation_and_azimuth(source, start_time, stop_time, coords).alt.value

        scan_list = []
        if source == "J2253+417":
            scan_list = source_J2253_417_scans
        elif  source == "EV_LAC":
            scan_list = source_EV_LAC_scans
        
        fig1, ax1 = plt.subplots(nrows=8, ncols=1, figsize=(16, 16), dpi=150, sharex=True)
        fig1.suptitle(source + " raw " + "dynamic spectrum")

        fig2, ax2 = plt.subplots(nrows=1, ncols=8, figsize=(16, 16), dpi=150, sharey=True)
        fig2.suptitle(source + " raw " + " IF plots ")

        fig3, ax3 = plt.subplots(nrows=8, ncols=1, figsize=(16, 16), dpi=150, sharex=True)
        fig3.suptitle(source + " banpass corrected " + " dynamic spectrum ")

        fig4, ax4 = plt.subplots(nrows=1, ncols=8, figsize=(16, 16), dpi=150, sharey=True)
        fig4.suptitle(source + " banpass corrected " + " IF plots ")

        fig5, ax5 = plt.subplots(nrows=1, ncols=8, figsize=(16, 16), dpi=150)
        fig5.suptitle(source + " raw " + " avg freq ")

        fig6, ax6 = plt.subplots(nrows=1, ncols=8, figsize=(16, 16), dpi=150)
        fig6.suptitle(source + " banpass corrected " + " avg freq ")

        fig7, ax7 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
        fig7.suptitle(source + " tsys ")

        fig8, ax8 = plt.subplots(nrows=8, ncols=1, figsize=(16, 16), dpi=150, sharex=True)
        fig8.suptitle(source + " gain corrected " + "dynamic spectrum")

        fig9, ax9 = plt.subplots(nrows=1, ncols=8, figsize=(16, 16), dpi=150, sharey=True)
        fig9.suptitle(source + " gain corrected " + " IF plots")
        
        fig10, ax10 = plt.subplots(nrows=1, ncols=8, figsize=(16, 16), dpi=150)
        fig10.suptitle(source + " gain corrected " + " avg freq")

        fig11, ax11 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150, sharex=True, sharey=True)
        fig11.suptitle(source + " model data")

        elevation_poly = [4.5642e-06, -0.0007491198, 1.0132377024]
        dfpu = [0.085300, 0.082000]

        for spw in range(0, nr_if):
            spw_table = taql('select from $source_table where DATA_DESC_ID == ' + str(spw))
            spw_data = spw_table.getcol("DATA")
            
            flag = 1-spw_table.getcol("FLAG")

            spw_data = np.ma.array(spw_data, mask=flag)
            
            ch_rr = tsys_ch_map[spw][0] - 1 
            ch_ll = tsys_ch_map[spw][1] - 1
            t_sys_rr = tsys[source][str(ch_rr)]
            t_sys_ll = tsys[source][str(ch_ll)]
            
            print(spw, "t_sys_rr", np.mean(t_sys_rr))
            print(spw, "t_sys_ll", np.mean(t_sys_ll))

            ax7.scatter(np.linspace(0, len(t_sys_rr), len(t_sys_rr)), t_sys_rr)
            ax7.scatter(np.linspace(0, len(t_sys_ll), len(t_sys_ll)), t_sys_ll)
            
            stokes_i_raw = np.abs(np.sum(spw_data[:,:,(0,3)],axis=-1))
            flag_raw = np.abs(np.sum(spw_data.mask[:,:,(0,3)],axis=-1))
            stokes_i_raw = np.ma.array(stokes_i_raw, mask=flag_raw)
            f_spw = frequencies_list_for_each_if[spw]
            

            im1 = ax1[spw].imshow(stokes_i_raw.T, aspect="auto",
                     extent=[np.min(time), np.max(time), np.min(f_spw), np.max(f_spw)],
                     vmin = np.percentile(stokes_i_raw.T, 1), vmax = np.percentile(stokes_i_raw.T, 95))

            ax2[spw].scatter(frequencies_list_for_each_if[spw], np.median(stokes_i_raw, axis=0))

            ll_corrected_banpass = []
            rr_corrected_banpass = []
            
            for scan in scan_list:
                scan_table = taql('select from $spw_table where SCAN_NUMBER == ' + str(scan))
                scan_data = spw_table.getcol("DATA")
                flag = 1-spw_table.getcol("FLAG")
                scan_data = np.ma.array(scan_data, mask=flag)
            
                ll = np.abs(scan_data[:,:,(0)])
                rr = np.abs(scan_data[:,:,(1)])
                banpass_correction_ll = np.median(ll[:1000, :], axis=0)
                banpass_correction_rr = np.median(ll[:1000, :], axis=0)

                ll_corrected_banpass.append(ll/banpass_correction_ll)
                rr_corrected_banpass.append(rr/banpass_correction_rr)
                                
            ll_corrected_banpass = np.vstack(ll_corrected_banpass)
            rr_corrected_banpass = np.vstack(rr_corrected_banpass)
            
            #ll_ = np.abs(scan_data[:,:,(0)])
            #rr_ = np.abs(scan_data[:,:,(1)])
            #ll_corrected_banpass_ = ll_/np.median(ll_[:1000, :], axis=0)
            #rr_corrected_banpass_ = rr_/np.median(rr_[:1000, :], axis=0)
            
            #print(rr_corrected_banpass_.shape, rr_corrected_banpass.shape)
            #sys.exit() 
            stokes_i_banpass_corrected = (ll_corrected_banpass+rr_corrected_banpass)/2
            nr_points = stokes_i_banpass_corrected.shape[0]
            
            elevation_intrep = np.interp(np.linspace(0, 1, nr_points),
                             np.linspace(0, 1, len(elevation)), elevation) 
            
            elevation_correction_left = []
            elevation_correction_right = []
            for el in elevation_intrep:
                elevation_correction_left.append((dfpu[0] * np.polyval(elevation_poly, el)))
                elevation_correction_right.append((dfpu[1] * np.polyval(elevation_poly, el)))
           
            elevation_correction_left = np.array(elevation_correction_left)
            elevation_correction_right = np.array(elevation_correction_right)
           
            t_sys_rr_polyfit = np.polyfit(np.linspace(0, 1, len(t_sys_rr)), t_sys_rr, 1) 
            t_sys_ll_polyfit = np.polyfit(np.linspace(0, 1, len(t_sys_ll)), t_sys_ll, 1)

            t_sys_rr_poly1d = np.poly1d(t_sys_rr_polyfit) 
            t_sys_ll_poly1d = np.poly1d(t_sys_ll_polyfit)

            t_sys_rr_polyval = t_sys_rr_poly1d(np.linspace(0, 1, nr_points))
            t_sys_ll_polyval = t_sys_ll_poly1d(np.linspace(0, 1, nr_points))
             
            amplitude_correction_ll = t_sys_ll_polyval  / elevation_correction_left
            amplitude_correction_rr = t_sys_rr_polyval / elevation_correction_right 

            ll_corrected_gain = ll_corrected_banpass * amplitude_correction_ll[:, None]
            rr_corrected_gain = rr_corrected_banpass * amplitude_correction_rr[:, None]       

            stokes_i_gain_corrected = ll_corrected_gain # (ll_corrected_gain+rr_corrected_gain)/2
            
            im2 = ax3[spw].imshow(stokes_i_banpass_corrected.T, aspect="auto",
                     extent=[np.min(time), np.max(time), np.min(f_spw), np.max(f_spw)],
                     vmin = np.percentile(stokes_i_banpass_corrected.T, 1), vmax = np.percentile(stokes_i_banpass_corrected.T, 95))
            
            ax4[spw].scatter(frequencies_list_for_each_if[spw], np.median(stokes_i_banpass_corrected, axis=0))
            ax4[spw].scatter(frequencies_list_for_each_if[spw], np.average(stokes_i_banpass_corrected, axis=0))
            
            avg_freq_raw =  np.median(stokes_i_raw, axis=1)
            
            ax5[spw].scatter(np.linspace(0,1, len(avg_freq_raw)), avg_freq_raw)

            avg_freq_banpass_corrected =  np.median(stokes_i_banpass_corrected, axis=1)
            avg_freq_banpass_corrected_a =  np.average(stokes_i_banpass_corrected, axis=1)
            ax6[spw].scatter(np.linspace(0,1, len(avg_freq_banpass_corrected)), avg_freq_banpass_corrected)
            ax6[spw].scatter(np.linspace(0,1, len(avg_freq_banpass_corrected)), avg_freq_banpass_corrected_a)

            im3 = ax8[spw].imshow(stokes_i_gain_corrected.T, aspect="auto",
                     extent=[np.min(time), np.max(time), np.min(f_spw), np.max(f_spw)],
                     vmin = np.percentile(stokes_i_gain_corrected.T, 1), vmax = np.percentile(stokes_i_gain_corrected.T, 95))
            
            ax9[spw].scatter(frequencies_list_for_each_if[spw], np.median(stokes_i_gain_corrected, axis=0))
            
            avg_freq_gain_corrected = np.median(stokes_i_gain_corrected, axis=1)
            ax10[spw].scatter(np.linspace(0,1, len(avg_freq_gain_corrected)), avg_freq_gain_corrected, label="data")

            if calibrator:
                polyfit = np.polyfit(np.linspace(0, 1, len(avg_freq_gain_corrected)), avg_freq_gain_corrected, 1) 
                poly1d = np.poly1d(polyfit) 
                output_calibrator_poly1d.append(poly1d)

                polyval = poly1d(np.linspace(0, 1, len(avg_freq_gain_corrected)))
                output = avg_freq_gain_corrected / polyval
                output_calibrator.append(output)
                ax11.scatter(np.linspace(0,1, len(output)), output, label=str(spw))
                ax10[spw].scatter(np.linspace(0,1, len(avg_freq_gain_corrected)), polyval, label="model")
                ax10[spw].scatter(np.linspace(0,1, len(avg_freq_gain_corrected)), output, label="result")

            else:
                calibrator_poly1d = output_calibrator_poly1d[spw]
                polyval = calibrator_poly1d(np.linspace(0, 1, len(avg_freq_gain_corrected)))

                output = avg_freq_gain_corrected / polyval
                output_target.append(output)
                ax11.scatter(np.linspace(0,1, len(output)), output, label=str(spw))
                ax10[spw].scatter(np.linspace(0,1, len(avg_freq_gain_corrected)), polyval, label="model")
                ax10[spw].scatter(np.linspace(0,1, len(avg_freq_gain_corrected)), output, label="result")
        
        ax11.legend()
        ax10[spw].legend()

        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        fig4.tight_layout()
        fig5.tight_layout()
        fig6.tight_layout()
        fig7.tight_layout()
        fig8.tight_layout()
        fig9.tight_layout()
        fig10.tight_layout()
        fig11.tight_layout()

        fig1.subplots_adjust(wspace=0, hspace=0)
        fig2.subplots_adjust(wspace=0, hspace=0)
        fig3.subplots_adjust(wspace=0, hspace=0)
        fig4.subplots_adjust(wspace=0, hspace=0)
        fig5.subplots_adjust(wspace=0, hspace=0)
        fig6.subplots_adjust(wspace=0, hspace=0)
        fig7.subplots_adjust(wspace=0, hspace=0)
        fig8.subplots_adjust(wspace=0, hspace=0)
        fig9.subplots_adjust(wspace=0, hspace=0)
        fig10.subplots_adjust(wspace=0, hspace=0)
        fig11.subplots_adjust(wspace=0, hspace=0)
  
        fig1.savefig("_".join([source, "raw", "dynamic_spectrum"]) + ".png")
        fig2.savefig("_".join([source, "raw", "IF_plots"]) + ".png" )
        fig3.savefig("_".join([source, "banpass_corrected", "dynamic_spectrum"]) + ".png")
        fig4.savefig("_".join([source, "banpass_corrected", "IF_plots"]) + ".png")
        fig5.savefig("_".join([source, "raw", "avg_freq"]) + ".png")
        fig6.savefig("_".join([source, "banpass_corrected", "avg_freq"]) + ".png")
        fig7.savefig("_".join([source, "tsys"]) + ".png")
        fig8.savefig("_".join([source, "gain_corrected", "dynamic_spectrum"]) + ".png")
        fig9.savefig("_".join([source, "gain_corrected", "IF_plots"]) + ".png")
        fig10.savefig("_".join([source, "gain_corrected", "avg_freq"]) + ".png")
        fig11.savefig("_".join([source, "model_data"]) + ".png")

        plt.close('all')

        if calibrator:
            calibrator = False


    fig12, ax12 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150, sharex=True, sharey=True)
    fig12.suptitle("Final results calibrator")
     
    fig13, ax13 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150, sharex=True, sharey=True)
    fig13.suptitle("Final results target")
    
    output_calibrator = np.array(output_calibrator)
    output_target = np.array(output_target)
    
    result_calibrator = np.zeros((len(output_calibrator[0])))
    result_target = np.zeros((len(output_target[0])))

    for i in range(0, nr_if):
        if i == 2 and i == 3:
            continue 

        result_calibrator += output_calibrator[i]/6
        result_target += output_target[i]/6

    ax12.scatter(np.linspace(0,1, len(result_calibrator)), result_calibrator)
    ax13.scatter(np.linspace(0,1, len(result_target)), result_target)
    
    fig12.savefig("Final_results_calibrator.png")
    fig13.savefig("Final_results_target.png")
    plt.close('all')           
    #plt.show()


if __name__ ==  "__main__":
    main()
    sys.exit(0)

