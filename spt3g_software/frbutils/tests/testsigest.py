#!/usr/bin/env python
from spt3g import core, dfmux, mapmaker, frbutils, todfilter
import numpy as np
import unittest, copy

def test_sig_est():
        ts_len = 101
        model_width = 21
        poly_order = 1

        inject_1_index = 33
        inject_2_index = 66

        ts_1 = core.G3Timestream(np.zeros( ts_len ))
        ts_2 = core.G3Timestream(np.zeros( ts_len ))

        ts_1[inject_1_index] = 10
        ts_1[inject_1_index+2] = 9
        ts_1[inject_1_index-2] = 9

        ts_2[inject_2_index] = 10

        ts_map = core.G3TimestreamMap()
        ts_map['ts_1'] = ts_1
        ts_map['ts_2'] = ts_2
        ts_map['ts_3'] = ts_1

        frame = core.G3Frame( core.G3FrameType.Scan )
        frame['TsMap'] = ts_map


        variance = core.G3MapDouble()
        variance['ts_1'] = 1
        variance['ts_2'] = 1
        variance['ts_3'] = 1
        frame['TsVariance'] = variance

        ts_key_to_use = 'TsMap'
        plf_base = frbutils.PolyLikelihoodFiller(
            model_len = model_width,
            poly_order = poly_order,
            include_heaviside = True,
            include_delta = False,
            ts_key = ts_key_to_use,
            variance_key = 'TsVariance',
            amp_map_output_key = 'FitDeltaAmplitude_baseline',
            hamp_map_output_key = 'FitHeaviAmplitude_baseline',#shouldn't be stored
            loglike_output_key = 'BaselineLogLike')
        
        plf_mod = frbutils.PolyLikelihoodFiller(
            model_len = model_width,
            poly_order = poly_order,
            include_heaviside = True,
            include_delta = True,
            ts_key = ts_key_to_use,
            variance_key = 'TsVariance',
            amp_map_output_key = 'FitDeltaAmplitude',
            hamp_map_output_key = 'FitHeaviAmplitude',
            loglike_output_key = 'ModelLogLike')

        ev_hunter = frbutils.DeltaEventHunter(
            ll_model_key = 'ModelLogLike',
            ll_base_key = 'BaselineLogLike',
            trigger_thresh = 10,
            other_det_thresh = 10,
            min_distance = 10,
            output_event_key = 'FrbEvents',
            fit_amp_key = 'FitDeltaAmplitude',
            fit_hamp_key = 'FitHeaviAmplitude',
            search_width = 1
        )

        plf_base(frame)
        plf_mod(frame)
        ev_hunter(frame)

        if False:
            import pylab as pl

            pl.clf()            
            pl.plot(frame['TsMap']['ts_1'])
            pl.title("Input Timestream")
            pl.xlabel("Sample")
            pl.savefig("/home/nlharr/Thesis/FrbPlots/LookinInputSignal.png")
            pl.clf()
            pl.title("Sample Significance")
            pl.xlabel("Sample")
            pl.plot(frame['ModelLogLike']['ts_1'] - frame['BaselineLogLike']['ts_1'])
            pl.savefig("/home/nlharr/Thesis/FrbPlots/LookinSignficance.png")

            #import pylab as pl; pl.plot(frame['FitDeltaAmplitude']['ts_1']); pl.plot(frame['FitHeaviAmplitude']['ts_1']); pl.show()
            #pl.ion()
            #import pdb, rlcompleter
            #pdb.Pdb.complete = rlcompleter.Completer(locals()).complete
            #pdb.set_trace()

        frb_evs = frame['FrbEvents']
        assert(len(frb_evs) == 2)
        assert(len(frb_evs[0].det_info)== 2)
        assert(len(frb_evs[1].det_info)== 1)

        assert(frb_evs[0].det_info[0].bid== 'ts_1')
        assert(frb_evs[0].det_info[1].bid == 'ts_3')
        assert(frb_evs[1].det_info[0].bid== 'ts_2')

        assert(frb_evs[0].scan_index== inject_1_index)
        assert(frb_evs[1].scan_index== inject_2_index)



def test_sig_inject():
        from spt3g.frbutils.impulseinjections import InjectFrbSignal
        from spt3g.frbutils import FastTransientSignalBeam 

        import os.path, pickle
        fname = '/home/nlharr/frb_side_products/fast_transient_signal_conversion.pkl'
        if not os.path.isfile(fname):
            print("one test file not found")
            return

        fts = pickle.load(open(fname))
        os.path.isfile(fname) 
        ts_len = 101
        model_width = 21
        poly_order = 1

        ts_1 = core.G3Timestream(np.zeros( ts_len ))
        ts_2 = core.G3Timestream(np.zeros( ts_len ))

        ts_map = core.G3TimestreamMap()
        ts_map['ts.X'] = ts_1
        ts_map['ts.Y'] = ts_2

        frame = core.G3Frame( core.G3FrameType.Scan )
        frame['TsMapEmpty'] = ts_map

        ts_key_to_use = 'TsMap'

        plf_base = frbutils.PolyLikelihoodFiller(
            model_len = model_width,
            poly_order = poly_order,
            include_heaviside = True,
            include_delta = False,
            ts_key = ts_key_to_use,
            variance_key = 'TsVariance',
            amp_map_output_key = 'FitDeltaAmplitude_baseline',
            hamp_map_output_key = 'FitHeaviAmplitude_baseline',#shouldn't be stored
            loglike_output_key = 'BaselineLogLike')
        
        plf_mod = frbutils.PolyLikelihoodFiller(
            model_len = model_width,
            poly_order = poly_order,
            include_heaviside = True,
            include_delta = True,
            ts_key = ts_key_to_use,
            variance_key = 'TsVariance',
            amp_map_output_key = 'FitDeltaAmplitude',
            hamp_map_output_key = 'FitHeaviAmplitude',
            loglike_output_key = 'ModelLogLike')

        ev_hunter = frbutils.DeltaEventHunter(
            ll_model_key = 'ModelLogLike',
            ll_base_key = 'BaselineLogLike',
            trigger_thresh = 10,
            other_det_thresh = 10,
            min_distance = 20,
            output_event_key = 'FrbEvents',
            fit_amp_key = 'FitDeltaAmplitude',
            fit_hamp_key = 'FitHeaviAmplitude',
            search_width = 1
        )

        InjectFrbSignal(frame, ts_key = 'TsMapEmpty',
                        out_ts_key = 'TsMap',
                        time_scale = 0.01, curve_type = 3,
                        fluence = 1e6, fts = fts, 
                        inject_index = 42
        )

        var_adder = todfilter.VarianceAdder(
            ts_key = 'TsMap',
            variance_output_key= 'TsVariance')
        
        var_adder(frame)
        plf_base(frame)
        plf_mod(frame)
        ev_hunter(frame)

        if 0:
            import pdb, rlcompleter
            pdb.Pdb.complete = rlcompleter.Completer(locals()).complete
            pdb.set_trace()

        frb_evs = frame['FrbEvents']
        assert(frb_evs[0].scan_index == frame['InjectedIndex'])

def test_individual_sig_est():
        ts_len = 101
        model_width = 21
        poly_order = 1

        inject_1_index = 33
        inject_2_index = 66

        ts_1 = core.G3Timestream(np.zeros( ts_len ))

        ts_1[inject_1_index] = 10
        ts_1[inject_1_index+2] = 9
        ts_1[inject_1_index-2] = 9


        ts_map = core.G3TimestreamMap()
        ts_map['ts_1'] = ts_1

        frame = core.G3Frame( core.G3FrameType.Scan )
        frame['TsMap'] = ts_map

        variance = core.G3MapDouble()
        variance['ts_1'] = 1
        frame['TsVariance'] = variance
        ts_key_to_use = 'TsMap'
        plf_mod = frbutils.PolyLikelihoodFiller(
            model_len = model_width,
            poly_order = poly_order,
            include_heaviside = True,
            include_delta = True,
            ts_key = ts_key_to_use,
            variance_key = 'TsVariance',
            amp_map_output_key = 'FitDeltaAmplitude',
            hamp_map_output_key = 'FitHeaviAmplitude',
            loglike_output_key = 'ModelLogLike')

        plf_mod(frame)

        amp_map = core.G3TimestreamMap()
        heavi_amp_map = core.G3TimestreamMap()

        for i in range(20, 80):
            ll_model_map = core.G3TimestreamMap()
            frbutils.get_ts_poly_delta_heavi_ll(ts_map, variance,
                                       model_width, poly_order,
                                       True, True, #heavi, delta,
                                       ll_model_map, amp_map, heavi_amp_map,
                                       i)
            assert( ll_model_map['ts_1'][0] == frame['ModelLogLike']['ts_1'][i] )

        
if __name__ == '__main__':
    test_sig_est()
    #test_sig_inject()
    test_individual_sig_est()
