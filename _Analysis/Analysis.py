import os
import mne
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import fnmatch



##############################################################################################
####                   Create Functions for vmrk Clean-Up                                 ####
##############################################################################################

def get_number_string(lyne):
    if line == '':
        number_str = "1"  # to get out of bottom while loop when file complete
    else:
        number_str = lyne.replace(' ', '').split(',')[1].strip().strip('S').strip('R')  # this takes the number
    return(number_str)


def do_vmrk_clean_up(participant):

    filedir = "/Volumes/LaCie/Experiments/MMIOR/Hannah/_Data/forMNE"

    filename = "multimodal_ior_follow_up_%s.vmrk" % participant

    os.chdir(filedir)  # change to file's folder

    # rename original it hasn't already been 'cleaned up'
    already_done = os.path.isfile("original_vmrks/original_" + filename)

    if not already_done:
        os.rename(filename, "original_vmrks/original_" + filename)

        newFile = open(filename, "w")  # open file in which error-free vmrk info will be put in
        oldFile = open("original_vmrks/original_" + filename, 'r')
        done = False
        dataStarted = False

        def get_number_string(lyne):
            if line == '':
                number_str = "1"  # to get out of bottom while loop when file complete
            else:
                number_str = lyne.replace(' ', '').split(',')[1].strip().strip('S').strip('R')  # this takes the number
            return (number_str)

        while not done:
            if not dataStarted:
                line = oldFile.readline()
                if line[0:3] == "Mk2":
                    dataStarted = True
                    number_string = get_number_string(line)
                else:
                    newFile.write(line)
            if dataStarted:
                if line == '':
                    done = True
                else:
                    if number_string == "1":
                        lines = line
                        line = oldFile.readline()
                        number_string = get_number_string(line)
                        while number_string != "1":
                            if number_string in ["99", "98", "97", "43", "42"]:
                                break
                            else:
                                lines += line
                            line = oldFile.readline()
                            number_string = get_number_string(line)
                        if number_string == "1":
                            newFile.write(lines)
                        else:
                            while number_string != "1":
                                line = oldFile.readline()
                                number_string = get_number_string(line)
                    else:
                        print("ERROR")

        newFile.close()
        oldFile.close()
    else:
        print("I HAVE ALREADY 'CLEANED UP' THIS FILE")



##############################################################################################
####                                Load Data                                             ####
##############################################################################################

# Participant
participants = ['e08', 'e10', 'e14']

dont_plot = True

for participant in participants:

    do_vmrk_clean_up(participant)

    os.chdir("/Volumes/LaCie/Experiments/MMIOR/Hannah/_Data")

    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if fnmatch.fnmatch(name, '*%s.vhdr'%participant):
                file = os.path.join(root, name)

    # for one participant
    raw = mne.io.read_raw_brainvision(file, preload=True)

    directory = '/Users/ghislaindentremont/Documents/Experiments/Multimodal_IOR/Hannah/P_analysis/%s'%participant
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists('%s/P_preprocessing'%directory):
        os.makedirs('%s/P_preprocessing'%directory)
    if not os.path.exists('%s/P_averages'%directory):
        os.makedirs('%s/P_averages'%directory)
    if not os.path.exists('%s/P_AR'%directory):
        os.makedirs('%s/P_AR'%directory)

    # remove Aux
    raw.drop_channels(['Aux1'])  # AUX1 = microsensor, o.w. nothing
    if  'Aux2' in raw.info['ch_names']:
        raw.drop_channels(['Aux2', 'Aux3', 'Aux4', 'Aux5', 'Aux6', 'Aux7', 'Aux8'])


    #--------------------------------- Look for bad channels ------------------------------------#
    raw.plot(n_channels = 64, scalings = dict(eeg = 200e-6), show = False)
    plt.savefig('%s/P_preprocessing/raw_waveforms_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()
    #--------------------------------- Look for bad channels ------------------------------------#


    #------------------------------------ Visualize Events --------------------------------------#
    events = mne.find_events(raw, stim_channel = 'STI 014', output = 'onset')
    mne.viz.plot_events(events, raw.info['sfreq'], raw.first_samp, show = False)
    plt.savefig('%s/P_preprocessing/events_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()
    #------------------------------------ Visualize Events --------------------------------------#


    #------------------------------------ AR on all channels ------------------------------------#
    picks_all = mne.pick_types(raw.info, meg=False, eeg=True, eog=False, exclude='bads')
    event_id_all = {
        'LV/LV': 20
        , 'LV/RV': 22
        , 'LT/LV': 24
        , 'LT/RV': 26
        , 'RV/LV': 28
        , 'RV/RV': 30
        , 'RT/LV': 32
        , 'RT/RV': 34

        , '700/LV/LV': 60
        , '700/LV/RV': 62
        , '700/LT/LV': 64
        , '700/LT/RV': 66
        , '700/RV/LV': 68
        , '700/RV/RV': 70
        , '700/RT/LV': 72
        , '700/RT/RV': 74
    }
    tmin, tmax = -0.2, 0.5
    epochs_params_all = dict(picks = picks_all, events=events, event_id=event_id_all, tmin=tmin, tmax=tmax)
    epochs_all = mne.Epochs(
        raw
        , **epochs_params_all
        , add_eeg_ref = False
        , baseline = (-0.1,0)
    )

    # artifact rejection
    epochs_all.drop_bad(reject = dict(eeg=100e-5), flat = dict(eeg = 100e-6))

    # percentage rejected by channel
    epochs_all.plot_drop_log(show = False)
    plt.savefig('%s/P_preprocessing/drop_log_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()
    #------------------------------------ AR on all channels ------------------------------------#


    #--------------------------------- Label Bad for each P -------------------------------------#
    if participant == "e01":
        raw.info['bads'] = []
    elif participant == "e02":
        raw.info['bads'] = []
    elif participant == "e03":
        raw.info['bads'] = ['Ch41']
    elif participant == "e04":
        raw.info['bads'] = []
    elif participant == "e05":
        raw.info['bads'] = ['Ch5']
    #--------------------------------- Label Bad for each P -------------------------------------#


    #--------------------------------------- See Montage ----------------------------------------#
    ch_names = [
        # greens (1-32)
        'Fp1'
        , 'Fz'
        , 'F3'
        , 'F7'
        , 'FT9'
        , 'FC5'
        , 'FC1'
        , 'C3'
        , 'T7'
        , 'TP9'
        , 'CP5'
        , 'CP1'
        , 'Pz'
        , 'P3'
        , 'P7'
        , 'O1'
        , 'Oz'
        , 'O2'
        , 'P4'
        , 'P8'
        , 'TP10'
        , 'CP6'
        , 'CP2'
        , 'Cz'
        , 'C4'
        , 'T8'
        , 'FT10'
        , 'FC6'
        , 'FC2'
        , 'F4'
        , 'F8'
        , 'Fp2'
        # yellows (33-64)
        , 'AF7'
        , 'AF3'
        , 'AFz'
        , 'F1'
        , 'F5'
        , 'FT7'
        , 'FC3'
        , 'FCz'
        , 'C1'
        , 'C5'
        , 'TP7'
        , 'CP3'
        , 'P1'
        , 'P5'
        , 'PO7'
        , 'PO3'
        , 'POz'
        , 'PO4'
        , 'PO8'
        , 'P6'
        , 'P2'
        , 'CPz'
        , 'CP4'
        , 'TP8'
        , 'C6'
        , 'C2'
        , 'FC4'
        , 'FT8'
        , 'F6'
        , 'F2'
        , 'AF4'
        , 'AF8'
        ]

    # get montage for interpolation and visualization
    montage = mne.channels.read_montage(kind ="easycap-M1", ch_names = ch_names)
    mne.viz.plot_montage(montage, show_names = True, show = False)
    if dont_plot:
        plt.close()
    #--------------------------------------- See Montage ----------------------------------------#


    #--------------------------------------- Set Montage ----------------------------------------#
    original_ch_names = raw.info['ch_names'][0:64]
    ch_names_dict = dict(zip(original_ch_names, ch_names))
    raw.rename_channels(mapping = ch_names_dict)

    raw.set_montage(montage)
    #--------------------------------------- Set Montage ----------------------------------------#


    #------------------------------------ Interpolate Bads --------------------------------------#
    if len(raw.info['bads']) > 0:
        raw.interpolate_bads()
        raw.plot(n_channels = 64, scalings = dict(eeg = 200e-6), show = False)
        plt.savefig('%s/P_preprocessing/waveforms_post_interpolation_%s.png'%(directory, participant) )
        if dont_plot:
            plt.close()
    #------------------------------------ Interpolate Bads --------------------------------------#





    ##############################################################################################
    ####                                Re-Reference                                          ####
    ##############################################################################################

    epochs_params_test = dict(picks = picks_all, events=events, event_id=event_id_all, tmin=tmin, tmax=tmax)


    #------------------------- Show Effect of Reference on Evoked -------------------------------#
    fig, ax = plt.subplots(2,1, sharex = True)

    # no reference
    raw_no_ref, _ = mne.io.set_eeg_reference(raw, [])
    evoked_no_ref = mne.Epochs(raw_no_ref, **epochs_params_test).average()
    del raw_no_ref  # save memory

    evoked_no_ref.plot(axes = ax[0], titles=dict(eeg='EEG Original reference'), show=False)

    # Average reference
    raw_ref, _ = mne.io.set_eeg_reference(raw)
    evoked_ref = mne.Epochs(raw_ref, **epochs_params_test).average()
    del raw_ref  # save memory

    evoked_ref.plot(axes = ax[1], titles=dict(eeg='EEG Average reference'), show=False)
    plt.savefig('%s/P_preprocessing/waveforms_pre_post_rereference_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()

    # # average of mastoid reference
    # raw_mast_ref, _ = mne.io.set_eeg_reference(raw, ['TP9', 'TP10'])
    # evoked_mast_ref = mne.Epochs(raw_mast_ref, **epochs_params_test).average()
    # del raw_mast_ref  # save memory
    #
    # evoked_mast_ref.plot(axes = ax[1], titles=dict(eeg='EEG Mastoid reference'))
    #------------------------- Show Effect of Reference on Evoked -------------------------------#


    #---------------------- Show Effect of Reference on Continuous ------------------------------#
    raw.plot(n_channels = 64, scalings = dict(eeg = 100e-6), show = False)
    plt.savefig('%s/P_preprocessing/waveforms_post_rereference_by_channel_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()
    #---------------------- Show Effect of Reference on Continuous ------------------------------#


    #------------------------ Apply and Rid Reference Projector ---------------------------------#
    raw.apply_proj()
    raw.info['projs'] = []
    #------------------------ Apply and Rid Reference Projector ---------------------------------#





    ##############################################################################################
    ####                                    Filters                                           ####
    ##############################################################################################

    #---------------------------------- Butterworth Filter --------------------------------------#
    # Pick a subset of channels
    picks = mne.pick_types(
        raw.info
        , meg = False
        , eeg = True
        , eog = False
        , selection =
        [ 'PO7', 'PO8']  # occipital(left/right)
        )

    # without filter
    raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (1,0,0), show = False)
    plt.savefig('%s/P_preprocessing/psd_pre_bandpass_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()

    # apply filter
    # applies zero-phase bandpass filter (default is 4th order butterworth, but using 2nd order here)
    iir_params = dict(order=2, ftype='butter')
    raw.filter(
        .1
        , 50
        , picks=picks
        , method = 'iir'
        , iir_params = iir_params
        , n_jobs = 4
    )
    raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (0,1,0), show = False)
    plt.savefig('%s/P_preprocessing/psd_post_bandpass_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()
    #---------------------------------- Butterworth Filter --------------------------------------#


    #---------------------------------- Notch Filter --------------------------------------------#
    # without filter
    raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (1,0,0), show = False)
    plt.savefig('%s/P_preprocessing/psd_pre_notch_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()

    # apply notch
    # zero phase notch filter
    raw.notch_filter(
        60
        , picks=picks
        , n_jobs = 4
        )
    raw.plot_psd(area_mode='range', tmax=10.0, picks=picks, color = (0,1,0), show = False)
    plt.savefig('%s/P_preprocessing/psd_post_notch_%s.png'%(directory, participant) )
    if dont_plot:
        plt.close()
    #---------------------------------- Notch Filter --------------------------------------------#





    ##############################################################################################
    ####                                   Functions                                          ####
    ##############################################################################################

    def get_evoked( raw, event_id, channel_name, tmin, tmax, reject_num, baseline = (-0.1,0), get_evoked = True, save_AR = False):
        reject = dict(eeg=reject_num)
        picks = mne.pick_types(
        raw.info
        , meg = False
        , eeg = True
        , eog = False
        , selection = channel_name
        )
        params = dict(
        picks =  picks
        , events=events  # global variable
        , event_id=event_id
        , tmin=tmin
        , tmax=tmax
        )
        epochs= mne.Epochs(
        raw
        , **params
        , add_eeg_ref = False
        , baseline = baseline
        )

        # artifact rejection
        epochs.drop_bad(reject = reject)

        if save_AR:
            # percentage rejected by channel
            if not dont_plot:
                epochs.plot_drop_log()
                plt.savefig('%s/P_AR/%s_%s_drop_log_%s.png' % (directory, channel_name[0], list(event_id)[0], participant))
            # # visualize by channel, by epoch
            # epochs.plot()
            f = open('%s/P_AR/%s_%s_AR_output_%s.txt' % (directory, channel_name[0], list(event_id)[0], participant), 'w')
            # get percentage of epochs dropped
            f.write("Total number of epochs: ")
            f.write(str(len(epochs)))
            f.write("\nPercentage of epochs dropped: ")
            f.write(str(epochs.drop_log_stats()))
            f.close()

        if get_evoked:
            evoked = epochs.average()
            # evoked.plot()
            # average channels
            sums = np.zeros(len(evoked.data[0]))
            for i in range(0, len(evoked.data)):
                sums = evoked.data[i] + sums

            avg = np.array( [sums/len(evoked.data)] )
            info = mne.create_info(
                ch_names = [' '.join(channel_name)]
                , sfreq = raw.info['sfreq']
                , ch_types = 'eeg'
                )
            avg_wave = mne.EvokedArray(avg, info, tmin=tmin)
            # avg_wave.plot()
            to_return = avg_wave
        else:
            to_return = epochs
        return to_return


    def concatenate_epochs( epoch1, epoch2, ch_name, reversal = False ):
        epoch1_arr = epoch1.get_data()
        epoch2_arr = epoch2.get_data()
        comb_arr = np.concatenate( [epoch1_arr, epoch2_arr] )

        comb_arr_sum = np.zeros(len(comb_arr[0][0]))
        for idx in range(0, len(comb_arr)):
            temp = comb_arr[idx][0]
            comb_arr_sum = comb_arr_sum + temp

        comb_arr_avg = np.array( [comb_arr_sum/len(comb_arr)] ) # get type 2d array

        if reversal:
            comb_arr_avg = comb_arr_avg * (-1)

        info_comb = mne.create_info(
            ch_names = [ch_name]
            , sfreq = epoch1.info['sfreq']
            , ch_types = 'eeg'
            )

        evoked = mne.EvokedArray(comb_arr_avg, info_comb, tmin=tmin)
        # evoked.plot()

        return evoked



    ##############################################################################################
    ####                          Visual Grand Average                                        ####
    ##############################################################################################

    vis_1000ms_id = {
        'LV/LV': 20
        , 'LV/RV': 22
        , 'LT/LV': 24
        , 'LT/RV': 26
        , 'RV/LV': 28
        , 'RV/RV': 30
        , 'RT/LV': 32
        , 'RT/RV': 34
    }

    evoked_grand_avg_vis_1000ms = get_evoked(raw, vis_1000ms_id, ['PO7', 'PO8'], tmin, tmax, reject_num=100e-6)
    # evoked_grand_avg_vis.plot()
    grand_avg_vis_1000ms = evoked_grand_avg_vis_1000ms.data[0]

    vis_700ms_id = {
        '700/LV/LV': 60
        , '700/LV/RV': 62
        , '700/LT/LV': 64
        , '700/LT/RV': 66
        , '700/RV/LV': 68
        , '700/RV/RV': 70
        , '700/RT/LV': 72
        , '700/RT/RV': 74
    }

    evoked_grand_avg_vis_700ms = get_evoked(raw, vis_700ms_id, ['PO7', 'PO8'], tmin, tmax, reject_num=100e-6)
    # evoked_grand_avg_vis.plot()
    grand_avg_vis_700ms = evoked_grand_avg_vis_700ms.data[0]


    # ---------------------------------------- Plot ---------------------------------------------#
    samps = np.shape(grand_avg_vis_1000ms)[0]
    X = np.linspace(-200, 500, samps)  # 701

    f, ax = plt.subplots(1, 2, sharex=True, sharey=True)

    f.add_subplot(111, frameon=False)

    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel("time (ms)")
    plt.ylabel("voltage (V)", labelpad=15)

    ax[0].plot(X, grand_avg_vis_1000ms)
    ax[0].axhline(y=0, color='black')
    ax[0].axvline(x=0, linestyle='dashed', color='black')
    ax[0].set_title('1000 ms CTOA')

    ax[1].plot(X, grand_avg_vis_700ms)
    ax[1].axhline(y=0, color='black')
    ax[1].axvline(x=0, linestyle='dashed', color='black')
    ax[1].set_title('700 ms CTOA')

    ax[1].set_ylim(ax[1].get_ylim()[::-1])
    ax[1].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.gca().invert_yaxis()

    plt.savefig('%s/P_averages/grand_averages_%s.png'%(directory, participant) )

    if dont_plot:
        plt.close()
    # ---------------------------------------- Plot ---------------------------------------------#





    ##############################################################################################
    ####                                    Save Dataframe                                    ####
    ##############################################################################################

    dataset = list(
        zip(
            X
            , grand_avg_vis_1000ms.tolist()
            , grand_avg_vis_700ms.tolist()
            )
        )

    df = pd.DataFrame(
        data = dataset
        , columns =
            [
            'time (ms)'
            , "ZZZ0"  # where Z means collapsed across
            , "ZZZ7"
            ]
     )

    df.to_csv( '%s/P_averages/grand_averages_%s.csv'%(directory, participant) )




    ##############################################################################################
    ####                                    By Condition                                      ####
    ##############################################################################################


    #################################### Tactile/Visual/1000 #######################################

    #---------------------------------- Cued/Contra/TV --------------------------------------------#
    epochs_contra_LTLV = get_evoked( raw, {'LT/LV': 24}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RTRV = get_evoked( raw, {'RT/RV': 34}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_contra_TV = concatenate_epochs(epochs_contra_LTLV, epochs_contra_RTRV, 'cued contra TV')
    cued_contra_TV = evoked_cued_contra_TV.data[0]
    #---------------------------------- Cued/Contra/TV --------------------------------------------#


    #---------------------------------- Uncued/Contra/TV ------------------------------------------#
    epochs_contra_LTRV = get_evoked( raw, {'LT/RV': 26}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RTLV = get_evoked( raw, {'RT/LV': 32}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_contra_TV = concatenate_epochs(epochs_contra_LTRV, epochs_contra_RTLV, 'uncued contra TV')
    uncued_contra_TV = evoked_uncued_contra_TV.data[0]
    #---------------------------------- Uncued/Contra/TV ------------------------------------------#


    #---------------------------------- Cued/Ipsi/TV ----------------------------------------------#
    epochs_ipsi_LTLV = get_evoked( raw, {'LT/LV': 24}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RTRV = get_evoked( raw, {'RT/RV': 34}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_ipsi_TV = concatenate_epochs(epochs_ipsi_LTLV, epochs_ipsi_RTRV, 'cued ipsi TV')
    cued_ipsi_TV = evoked_cued_ipsi_TV.data[0]
    #---------------------------------- Cued/Ipsi/TV ----------------------------------------------#


    #---------------------------------- Uncued/Ipsi/TV --------------------------------------------#
    epochs_ipsi_LTRV = get_evoked( raw, {'LT/RV': 26}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RTLV = get_evoked( raw, {'RT/LV': 32}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_ipsi_TV = concatenate_epochs(epochs_ipsi_LTRV, epochs_ipsi_RTLV, 'uncued ipsi TV')
    uncued_ipsi_TV = evoked_uncued_ipsi_TV.data[0]
    #---------------------------------- Uncued/Ipsi/TV --------------------------------------------#




    #################################### Tactile/Visual/700 ########################################

    #---------------------------------- Cued/Contra/TV/700 ----------------------------------------#
    epochs_contra_LTLV_700 = get_evoked( raw, {'700/LT/LV': 64}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RTRV_700 = get_evoked( raw, {'700/RT/RV': 74}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_contra_TV_700 = concatenate_epochs(epochs_contra_LTLV_700, epochs_contra_RTRV_700, 'cued contra TV 700')
    cued_contra_TV_700 = evoked_cued_contra_TV_700.data[0]
    #---------------------------------- Cued/Contra/TV/700 ----------------------------------------#


    #---------------------------------- Uncued/Contra/TV/700 --------------------------------------#
    epochs_contra_LTRV_700 = get_evoked( raw, {'700/LT/RV': 66}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RTLV_700 = get_evoked( raw, {'700/RT/LV': 72}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_contra_TV_700 = concatenate_epochs(epochs_contra_LTRV_700, epochs_contra_RTLV_700, 'uncued contra TV 700')
    uncued_contra_TV_700 = evoked_uncued_contra_TV_700.data[0]
    #---------------------------------- Uncued/Contra/TV/700 --------------------------------------#


    #---------------------------------- Cued/Ipsi/TV/700 ------------------------------------------#
    epochs_ipsi_LTLV_700 = get_evoked( raw, {'700/LT/LV': 64}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RTRV_700 = get_evoked( raw, {'700/RT/RV': 74}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_ipsi_TV_700 = concatenate_epochs(epochs_ipsi_LTLV_700, epochs_ipsi_RTRV_700, 'cued ipsi TV 700')
    cued_ipsi_TV_700 = evoked_cued_ipsi_TV_700.data[0]
    #---------------------------------- Cued/Ipsi/TV/700 ------------------------------------------#


    #---------------------------------- Uncued/Ipsi/TV/700 ----------------------------------------#
    epochs_ipsi_LTRV_700 = get_evoked( raw, {'700/LT/RV': 66}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RTLV_700 = get_evoked( raw, {'700/RT/LV': 72}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_ipsi_TV_700 = concatenate_epochs(epochs_ipsi_LTRV_700, epochs_ipsi_RTLV_700, 'uncued ipsi TV 700')
    uncued_ipsi_TV_700 = evoked_uncued_ipsi_TV_700.data[0]
    #---------------------------------- Uncued/Ipsi/TV/700 ----------------------------------------#




    ####################################### Visual/Visual/1000 #####################################

    #---------------------------------- Cued/Contra/VV --------------------------------------------#
    epochs_contra_LVLV = get_evoked( raw, {'LV/LV': 20}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RVRV = get_evoked( raw, {'RV/RV': 30}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_contra_VV = concatenate_epochs(epochs_contra_LVLV, epochs_contra_RVRV, 'cued contra VV')
    cued_contra_VV = evoked_cued_contra_VV.data[0]
    #---------------------------------- Cued/Contra/VV --------------------------------------------#


    #---------------------------------- Uncued/Contra/VV ------------------------------------------#
    epochs_contra_LVRV = get_evoked( raw, {'LV/RV': 22}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RVLV = get_evoked( raw, {'RV/LV': 28}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_contra_VV = concatenate_epochs(epochs_contra_LVRV, epochs_contra_RVLV, 'uncued contra VV')
    uncued_contra_VV = evoked_uncued_contra_VV.data[0]
    #---------------------------------- Uncued/Contra/VV ------------------------------------------#


    #---------------------------------- Cued/Ipsi/VV ----------------------------------------------#
    epochs_ipsi_LVLV = get_evoked( raw, {'LV/LV': 20}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RVRV = get_evoked( raw, {'RV/RV': 30}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_ipsi_VV = concatenate_epochs(epochs_ipsi_LVLV, epochs_ipsi_RVRV, 'cued ipsi VV')
    cued_ipsi_VV = evoked_cued_ipsi_VV.data[0]
    #---------------------------------- Cued/Ipsi/VV ----------------------------------------------#


    #---------------------------------- Uncued/Ipsi/VV --------------------------------------------#
    epochs_ipsi_LVRV = get_evoked( raw, {'LV/RV': 22}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RVLV = get_evoked( raw, {'RV/LV': 28}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_ipsi_VV = concatenate_epochs(epochs_ipsi_LVRV, epochs_ipsi_RVLV, 'uncued ipsi VV')
    uncued_ipsi_VV = evoked_uncued_ipsi_VV.data[0]
    #---------------------------------- Uncued/Ipsi/VV --------------------------------------------#




    ####################################### Visual/Visual/700 ######################################

    #---------------------------------- Cued/Contra/VV/700 ----------------------------------------#
    epochs_contra_LVLV_700 = get_evoked( raw, {'700/LV/LV': 60}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RVRV_700 = get_evoked( raw, {'700/RV/RV': 70}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_contra_VV_700 = concatenate_epochs(epochs_contra_LVLV_700, epochs_contra_RVRV_700, 'cued contra VV 700')
    cued_contra_VV_700 = evoked_cued_contra_VV_700.data[0]
    #---------------------------------- Cued/Contra/VV/700 ----------------------------------------#


    #---------------------------------- Uncued/Contra/VV/700 --------------------------------------#
    epochs_contra_LVRV_700 = get_evoked( raw, {'700/LV/RV': 62}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_contra_RVLV_700 = get_evoked( raw, {'700/RV/LV': 68}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_contra_VV_700 = concatenate_epochs(epochs_contra_LVRV_700, epochs_contra_RVLV_700, 'uncued contra VV 700')
    uncued_contra_VV_700 = evoked_uncued_contra_VV_700.data[0]
    #---------------------------------- Uncued/Contra/VV/700 --------------------------------------#


    #---------------------------------- Cued/Ipsi/VV/700 ------------------------------------------#
    epochs_ipsi_LVLV_700 = get_evoked( raw, {'700/LV/LV': 60}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RVRV_700 = get_evoked( raw, {'700/RV/RV': 70}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_cued_ipsi_VV_700 = concatenate_epochs(epochs_ipsi_LVLV_700, epochs_ipsi_RVRV_700, 'cued ipsi VV 700')
    cued_ipsi_VV_700 = evoked_cued_ipsi_VV_700.data[0]
    #---------------------------------- Cued/Ipsi/VV/700 ------------------------------------------#


    #---------------------------------- Uncued/Ipsi/VV/700 ----------------------------------------#
    epochs_ipsi_LVRV_700 = get_evoked( raw, {'700/LV/RV': 62}, ['PO8'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    epochs_ipsi_RVLV_700 = get_evoked( raw, {'700/RV/LV': 68}, ['PO7'], tmin, tmax, reject_num = 100e-6, get_evoked = False )
    evoked_uncued_ipsi_VV_700 = concatenate_epochs(epochs_ipsi_LVRV_700, epochs_ipsi_RVLV_700, 'uncued ipsi VV 700')
    uncued_ipsi_VV_700 = evoked_uncued_ipsi_VV_700.data[0]
    #---------------------------------- Uncued/Ipsi/VV/700 ----------------------------------------#





    #------------------------------------ Plot Together -----------------------------------------#
    X = np.linspace(-200, 500, samps)

    f, ax = plt.subplots(2,4, sharex = True, sharey = True)

    ax111 = f.add_subplot(111, frameon = False)
    ax111.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax111.set_xlabel("time (ms)")
    ax111.set_ylabel("voltage (V)", labelpad = 15)

    # cue modality labels
    ax211 = f.add_subplot(211, frameon = False)
    ax211.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax211.set_ylabel("tactile cue", rotation = 270, labelpad = 20)
    ax211.yaxis.set_label_position("right")

    ax212 = f.add_subplot(212, frameon = False)
    ax212.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax212.set_ylabel("visual cue", rotation = 270, labelpad = 20)
    ax212.yaxis.set_label_position("right")

    # target modality labels
    ax121 = f.add_subplot(121, frameon = False)
    ax121.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax121.set_title("700 ms", y = 1.05)

    ax122 = f.add_subplot(122, frameon = False)
    ax122.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax122.set_title("1000 ms", y = 1.05)

    # contralaterality labels
    ax141 = f.add_subplot(141, frameon = False)
    ax141.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax141.set_title("contralateral")

    ax142 = f.add_subplot(142, frameon = False)
    ax142.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax142.set_title("ipsilateral")

    ax143 = f.add_subplot(143, frameon = False)
    ax143.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax143.set_title("contralateral")

    ax144 = f.add_subplot(144, frameon = False)
    ax144.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    ax144.set_title("ipsilateral")

    ax[0,0].plot(X, cued_contra_TV_700, label = 'cued')
    ax[0,0].plot(X, uncued_contra_TV_700, label = 'uncued')
    ax[0,0].axhline(y=0, color = 'black')
    ax[0,0].axvline(x=0, linestyle='dashed', color = 'black')
    ax[0,0].legend(prop = {'size':12})

    ax[0,1].plot(X, cued_ipsi_TV_700, label = 'cued')
    ax[0,1].plot(X, uncued_ipsi_TV_700, label = 'uncued')
    ax[0,1].axhline(y=0, color = 'black')
    ax[0,1].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,0].plot(X, cued_contra_VV_700, label = 'cued')
    ax[1,0].plot(X, uncued_contra_VV_700, label = 'uncued')
    ax[1,0].axhline(y=0, color = 'black')
    ax[1,0].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,1].plot(X, cued_ipsi_VV_700, label = 'cued')
    ax[1,1].plot(X, uncued_ipsi_VV_700, label = 'uncued')
    ax[1,1].axhline(y=0, color = 'black')
    ax[1,1].axvline(x=0, linestyle='dashed', color = 'black')

    ax[0,2].plot(X, cued_contra_TV, label = 'cued')
    ax[0,2].plot(X, uncued_contra_TV, label = 'uncued')
    ax[0,2].axhline(y=0, color = 'black')
    ax[0,2].axvline(x=0, linestyle='dashed', color = 'black')

    ax[0,3].plot(X, cued_ipsi_TV, label = 'cued')
    ax[0,3].plot(X, uncued_ipsi_TV, label = 'uncued')
    ax[0,3].axhline(y=0, color = 'black')
    ax[0,3].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,2].plot(X, cued_contra_VV, label = 'cued')
    ax[1,2].plot(X, uncued_contra_VV, label = 'uncued')
    ax[1,2].axhline(y=0, color = 'black')
    ax[1,2].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,3].plot(X, cued_ipsi_VV, label = 'cued')
    ax[1,3].plot(X, uncued_ipsi_VV, label = 'uncued')
    ax[1,3].axhline(y=0, color = 'black')
    ax[1,3].axvline(x=0, linestyle='dashed', color = 'black')

    ax[1,3].set_ylim(ax[1,3].get_ylim()[::-1])
    ax[1,3].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    f.set_size_inches(18.5, 10.5)

    plt.savefig( '%s/P_averages/condition_averages_%s.png'%(directory, participant) )

    if dont_plot:
        plt.close()
    #------------------------------------ Plot Together -----------------------------------------#



    ##############################################################################################
    ####                                    Save Dataframe                                    ####
    ##############################################################################################

    dataset = list(
        zip(
            X
            , cued_contra_TV_700.tolist()
            , uncued_contra_TV_700.tolist()
            , cued_ipsi_TV_700.tolist()
            , uncued_ipsi_TV_700.tolist()

            , cued_contra_VV_700.tolist()
            , uncued_contra_VV_700.tolist()
            , cued_ipsi_VV_700.tolist()
            , uncued_ipsi_VV_700.tolist()

            , cued_contra_TV.tolist()
            , uncued_contra_TV.tolist()
            , cued_ipsi_TV.tolist()
            , uncued_ipsi_TV.tolist()

            , cued_contra_VV.tolist()
            , uncued_contra_VV.tolist()
            , cued_ipsi_VV.tolist()
            , uncued_ipsi_VV.tolist()
            )
        )

    df = pd.DataFrame(
        data = dataset
        , columns =
            [
            'time (ms)'
            ,'CCT7'
            , 'UCT7'
            , 'CIT7'
            , 'UIT7'

            , 'CCV7'
            , 'UCV7'
            , 'CIV7'
            , 'UIV7'

            , 'CCT0'
            , 'UCT0'
            , 'CIT0'
            , 'UIT0'

            , 'CCV0'
            , 'UCV0'
            , 'CIV0'
            , 'UIV0'
            ]
     )

    df.to_csv( '%s/P_averages/condition_averages_%s.csv'%(directory, participant) )
