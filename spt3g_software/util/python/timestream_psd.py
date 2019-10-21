"""
This module contains functionality for plotting the for plotting many
detector PSDs together
"""
from spt3g.core import G3Units
from . import stats
import numpy as np
from itertools import islice
import copy
import astropy.stats as astrostats

def timestream_psd(timestream):
    from matplotlib.mlab import psd, window_hanning

    pad_to = 2 << int(np.log2(len(timestream)))
    sample_rate = timestream.sample_rate
    if sample_rate == float('inf'):
        raise RuntimeError("Sample rate is infinity")
    return np.array(psd(timestream, NFFT=len(timestream), pad_to=pad_to,
                        Fs=sample_rate, window=window_hanning))

def bolo_psd_plot(tsmap, low_f=None, high_f=None, bolo_order=None, vmin=None,
                  vmax=None):
    """
    Plot the PSD for each detector or given detectors over given frequenct range.
    """
    from matplotlib.mlab import psd, window_hanning

    import matplotlib.pyplot as plt
    if bolo_order is not None:
        subset = [(bolo, tsmap[bolo]) for bolo in bolo_order]
        fig, axes = plt.subplots(ncols=len(subset))
        for idx, (bolo, ts) in enumerate(subset):
            psd, freqs = timestream_psd(ts)
            if low_f is not None:
                start_idx = next(i for i in range(len(freqs)) if freqs[i]>low_f)
            else:
                start_idx = 0
            if high_f is not None and high_f<=freqs[-1]:
                stop_idx = next(i for i in range(start_idx, len(freqs)) if freqs[i]>high_f)
            else:
                stop_idx = len(freqs)-1
            ymin = freqs[start_idx]/G3Units.Hz
            ymax = freqs[stop_idx]/G3Units.Hz
            # Format psd for imshow
            psd = psd[start_idx:stop_idx][::-1]
            psd = psd.reshape(len(psd), 1)
            axes[idx].imshow(psd, aspect='auto', cmap=plt.get_cmap('jet'),
                             extent=[0, 1, ymin, ymax], vmin=vmin, vmax=vmax)
            axes[idx].get_xaxis().set_ticks([])
            if idx != 0:
                axes[idx].get_yaxis().set_ticks([])
            axes[idx].set_xlabel(bolo).set_rotation(70)
        axes[0].set_ylabel("Frequency (Hz)")
    else:
        fig, axes = plt.subplots()
        # Make bolos and psds have same order for clicky callback magic.
        bolos = tsmap.keys()
        freqs = timestream_psd(tsmap[bolos[0]])[1]
        if low_f is not None:
            start_idx = next(i for i in range(len(freqs)) if freqs[i]>low_f)
        else:
            start_idx = 0
        if high_f is not None and high_f<=freqs[-1]:
            stop_idx = next(i for i in range(start_idx, len(freqs)) if freqs[i]>high_f)
        else:
            stop_idx = len(freqs)-1
        psds = np.array([timestream_psd(tsmap[bolo])[0][start_idx:stop_idx]
                    for bolo in bolos])
        freqs = freqs[start_idx:stop_idx]

        tmp_psd = copy.deepcopy(psds).flatten()
        psd_std = astrostats.mad_std(list(tmp_psd))
        psd_mean = np.mean(tmp_psd)
        good_psds = [psd for psd in tmp_psd if
                (psd > psd_mean-psd_std and psd < psd_mean+psd_std)]
        psd_mean = np.mean(good_psds)
        if vmin is None:
            vmin = psd_mean-15*psd_std
        if vmax is None:
            vmax = psd_mean+15*psd_std

        heatmap = axes.pcolormesh(range(len(bolos)), freqs/G3Units.Hz,
                                  psds.transpose(), vmin=vmin, vmax=vmax)
        axes.set_xlim([0,len(bolos)])
        axes.set_ylim([freqs[0]/G3Units.Hz,freqs[-1]/G3Units.Hz])
        axes.set_ylabel("Frequency (Hz)")
        axes.set_xlabel("Bolometer")
        cb = fig.colorbar(heatmap)
        if tsmap[bolos[0]].units is not None:
            cb.set_label(r"PSD (${0}^2/Hz$)".format(tsmap[bolos[0]].units))

        text_props = {
            'mouse_loc': None,
            'text': None
            }
        def mouse_down(event):
            text_props['mouse_loc'] = (event.x, event.y)
        def mouse_up(event):
            if text_props['mouse_loc'] == (event.x, event.y):
                if text_props['text'] is None:
                    bolo = str(bolos[int(event.xdata)])
                    f = plt.gcf()
                    xdim, ydim = f.get_size_inches()*f.dpi
                    text_props['text'] = fig.text(event.x/xdim, event.y/ydim,
                                                  bolo, bbox={'alpha': 0.9})
                else:
                    text_props['text'].remove()
                    text_props['text'] = None
                event.canvas.draw()

        fig.canvas.mpl_connect('button_press_event', mouse_down)
        fig.canvas.mpl_connect('button_release_event', mouse_up)
        return fig
