print("""\n\n\n[USAGE]
Script designed to be used immediately after data collection from terminal
to quickly visualize waveforms without gaining information about experimental effects (potential bias).
No interpolation is performed and error trials are not removed\n\n\n""")

import os
import mne
import numpy as np
from matplotlib import pyplot as plt



##############################################################################################
####                                Load Data                                             ####
##############################################################################################

# get dir
filedir = raw_input("Where is the .eeg file found?\n>>> ")
# filedir = input()
# filedir = "/Users/ray/Experiments/Ghis-Multimodal-IOR/_EEG"

# Participant
participant = raw_input("What is the participants id (e.g. e01)?\n>>> ")
# participant = "e32"

# change dir
os.chdir(filedir)

# for one participant
raw = mne.io.read_raw_brainvision('multimodal_ior_follow_up_%s.vhdr' % participant, preload = True)

# remove Aux
raw.drop_channels(['Aux1']) # AUX = microsensor
if  'Aux2' in raw.info['ch_names']:
    raw.drop_channels(['Aux2', 'Aux3', 'Aux4', 'Aux5', 'Aux6', 'Aux7', 'Aux8'])

events = mne.find_events(raw, stim_channel = 'STI 014', output = 'onset')
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


#--------------------------------------- Set Montage ----------------------------------------#
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

original_ch_names = raw.info['ch_names'][0:64]
ch_names_dict = dict(zip(original_ch_names, ch_names))
raw.rename_channels(mapping = ch_names_dict)

raw.set_montage(montage)
#--------------------------------------- Set Montage ----------------------------------------#



##############################################################################################
####                                Re-Reference                                          ####
##############################################################################################

epochs_params_test = dict(picks = picks_all, events=events, event_id=event_id_all, tmin=tmin, tmax=tmax)


#------------------------- Show Effect of Reference on Evoked -------------------------------#
# no reference
mne.io.set_eeg_reference(raw, [])

# Average reference
mne.io.set_eeg_reference(raw)
#------------------------- Show Effect of Reference on Evoked -------------------------------#


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
    , selection = ['PO7', 'PO8']  # occipital(left/right)
    )

# apply filter
iir_params = dict(order=2, ftype='butter')
raw.filter(
    .1
    , 50
    , picks=picks
    , method = 'iir'
    , iir_params = iir_params
    # , n_jobs = 4
)
#---------------------------------- Butterworth Filter --------------------------------------#


#---------------------------------- Notch Filter --------------------------------------------#
# apply notch
raw.notch_filter(
    60  # get 60 
    , picks=picks
    # , n_jobs = 4
    )
#---------------------------------- Notch Filter --------------------------------------------#




##############################################################################################
####                                   Functions                                          ####
##############################################################################################

def get_evoked(raw, event_id, channel_name, tmin, tmax, reject_num, baseline=(-0.1, 0), get_evoked=True):
    reject = dict(eeg=reject_num)
    picks = mne.pick_types(
        raw.info
        , meg=False
        , eeg=True
        , eog=False
        , selection=channel_name
    )
    # params = dict(
    #     picks=picks
    #     , events=events  # global variable
    #     , event_id=event_id
    #     , tmin=tmin
    #     , tmax=tmax
    # )
    epochs = mne.Epochs(
        raw
        , picks=picks
        , events=events  # global variable
        , event_id=event_id
        , tmin=tmin
        , tmax=tmax
        , add_eeg_ref=False
        , baseline=baseline
    )
    # artifact rejection
    epochs.drop_bad(reject=reject)
    # # percentage rejected by channel
    # epochs.plot_drop_log()
    # # visualize by channel, by epoch
    # epochs.plot()
    if get_evoked:
        evoked = epochs.average()
        # evoked.plot()
        # average channels
        sums = np.zeros(len(evoked.data[0]))
        for i in range(0, len(evoked.data)):
            sums = evoked.data[i] + sums
        avg = np.array([sums / len(evoked.data)])
        info = mne.create_info(
            ch_names=[' '.join(channel_name)]
            , sfreq=raw.info['sfreq']
            , ch_types='eeg'
        )
        avg_wave = mne.EvokedArray(avg, info, tmin=tmin)
        # avg_wave.plot()
        to_return = avg_wave
    else:
        to_return = epochs
    return to_return
    
    


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

evoked_grand_avg_vis_1000ms = get_evoked( raw, vis_1000ms_id, ['PO7','PO8'], tmin, tmax, reject_num = 100e-6 )
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

evoked_grand_avg_vis_700ms = get_evoked( raw, vis_700ms_id, ['PO7','PO8'], tmin, tmax, reject_num = 100e-6 )
# evoked_grand_avg_vis.plot()
grand_avg_vis_700ms = evoked_grand_avg_vis_700ms.data[0]


#---------------------------------------- Plot ---------------------------------------------#
samps = np.shape(grand_avg_vis_1000ms)[0]
X = np.linspace(-200, 500, samps) # 701

f, ax = plt.subplots(1,2, sharex = True, sharey = True)

f.add_subplot(111, frameon = False)

# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel("time (ms)")
plt.ylabel("voltage (V)", labelpad = 15)

ax[0].plot(X, grand_avg_vis_1000ms)
ax[0].axhline(y=0, color = 'black')
ax[0].axvline(x=0, linestyle='dashed', color = 'black')
ax[0].set_title('1000 ms CTOA')

ax[1].plot(X, grand_avg_vis_700ms)
ax[1].axhline(y=0, color = 'black')
ax[1].axvline(x=0, linestyle='dashed', color = 'black')
ax[1].set_title('700 ms CTOA')

ax[1].set_ylim(ax[1].get_ylim()[::-1])
ax[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.gca().invert_yaxis()

plt.show()
#---------------------------------------- Plot ---------------------------------------------#

