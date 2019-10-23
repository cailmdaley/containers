# *************************************************************************** #
#
#   Main contents of this script:
#       1) A pipeline segment that combines three scan frames together,
#          applies certain user-supplied processing step
#          to the combined frame, and split that frame
#          back into three frames.
#       2) A pipeline that tests this pipeline segment.
#
# *************************************************************************** #



# =============================================================================
# Import modules
# -----------------------------------------------------------------------------

from spt3g import core
from spt3g import dfmux
from spt3g import gcp
from spt3g import coordinateutils
from spt3g import todfilter
import copy
import numpy as np


# =============================================================================





# =============================================================================
# Define pipeline modules and segments
# -----------------------------------------------------------------------------


# Things related to concatenation and deconcatenation
# -----------------------------------------------------------------------------

class ConcatenateScanFrames(object):
    
    '''
    This pipeline module concatenates the data stored in
    three consecutive scan frames, which correspond to a turnaround,
    a constant speed scan, and another turnaround, and returns one scan frame
    that has the concatenated data.
    
    '''
    
    def __init__(self):
        
        self.cache_of_frames        = []
        self.need_to_clear_cache    = False
        self.scan_frame_number      = 0
        self.cache_of_frame_numbers = []
    
    
    def copy_cartain_key_and_value_from_all_frames(
            self, individual_frames, key, combined_frame):
        
        for i in range(len(individual_frames)):
            new_key = "("+key+"FromScanFrame"+str(i)+"InTheCache"+")"
            combined_frame[new_key] = individual_frames[i][key]
    
    
    def concatenate_timestreams(self, frames, key):
        
        tss      = [frame[key] for frame in frames]
        long_ts  = core.G3Timestream.concatenate(tss)
        chk_lens = [ts.n_samples for ts in tss]
        # print("# Length of each chunk:")
        # print(chk_lens)
        return long_ts, chk_lens
    
    
    def concatenate_timestream_maps(self, frames, key, allowed_separation=0.6):
        
        tsms     = [frame[key] for frame in frames]
        chk_lens = [tsm.n_samples for tsm in tsms]
        long_tsm = core.G3TimestreamMap.concatenate(
                        tsms, ts_rounding_error=allowed_separation)
        return long_tsm, chk_lens
    
    
    def concatenate_g3_vectors(self, frames, key):
        
        vecs     = [frame[key] for frame in frames]
        vec_lens = [len(vec) for vec in vecs]
        long_vec = copy.copy(vecs[0])
        for vec in vecs[1:]:
            long_vec.extend(vec)
        return long_vec, vec_lens
    
    
    def concatenate_g3_map_vectors(self, frames, key):
        
        mvs = [frame[key] for frame in frames]
        mv  = type(frames[0][key])()
        mvl = core.G3MapVectorInt()
        for k in mvs[0].keys():
            vctrs     = [mv[k] for mv in mvs]
            vctr_lens = [len(vec) for vec in vctrs]
            mvl[k]    = vctr_lens
            long_vctr = copy.copy(vctrs[0])
            for vctr in vctrs[1:]:
                long_vctr.extend(vctr)
            mv[k] = long_vctr
        return mv, mvl
    
    
    def concatenate_frames(self, frames):
        
        all_keys = set()
        for frame in frames:
            all_keys = all_keys | set(frame.keys())
        
        combined_frame = core.G3Frame(core.G3FrameType.Scan)
        
        for key in all_keys:
            if key == "Turnaround":
                pass
            
            elif key in ["GCPFeatureBits", "ObservationID", "SourceName"]:
                # These keys store the same information in all scan frames,
                # so, the information from any frame can be copied.
                combined_frame[key] = frames[1][key]
            
            elif key in ["DfMuxHousekeeping", "OnlinePointingModel"]:
                # These keys also seem to store the same information in
                # all scan frames, but confirmation is needed.
                combined_frame[key] = frames[1][key]
                self.copy_cartain_key_and_value_from_all_frames(
                    frames, key, combined_frame)
            
            elif key in ["TrackerStatus", "TrackerPointing"]:
                # Concatenation can be easily done with "+" operator,
                # but deconcatenation seems a little more involved.
                # These keys seem to store read-only type of information,
                # so the information from each frame will be simply
                # copied in the combined frame and copied back later.
                concatenated_data = type(frames[1][key])()
                for frame in frames:
                    concatenated_data += frame[key]
                combined_frame[key] = concatenated_data
                self.copy_cartain_key_and_value_from_all_frames(
                    frames, key, combined_frame)
            
            elif isinstance(frame[key], core.G3Timestream):
                long_ts, chunk_lengths = \
                    self.concatenate_timestreams(frames, key)
                combined_frame[key] = long_ts
                combined_frame["("+key+"ChunkLensUsedForConcatenation"+")"] = \
                    core.G3VectorInt(chunk_lengths)
            
            elif isinstance(frame[key], core.G3TimestreamMap):
                if key == "BenchPosition":
                    sep = 1.0
                else:
                    sep = 0.6
                long_tsm, chunk_lengths = \
                    self.concatenate_timestream_maps(frames, key,
                                                     allowed_separation=sep)
                combined_frame[key] = long_tsm
                combined_frame["("+key+"ChunkLensUsedForConcatenation"+")"] = \
                    core.G3VectorInt(chunk_lengths)
            
            elif ("Vector"  in str(type(frame[key]))) and \
                 ("Map" not in str(type(frame[key]))):
                long_vec, vec_lengths = \
                    self.concatenate_g3_vectors(frames, key)
                combined_frame[key] = long_vec
                combined_frame["("+key+"ChunkLensUsedForConcatenation"+")"] = \
                    core.G3VectorInt(vec_lengths)
            
            elif ("Vector"  in str(type(frame[key]))) and \
                 ("Map"     in str(type(frame[key]))):
                long_mv, mv_lengths = \
                    self.concatenate_g3_map_vectors(frames, key)
                combined_frame[key] = long_mv
                combined_frame["("+key+"ChunkLensUsedForConcatenation"+")"] = \
                    core.G3VectorInt(mv_lengths)
            
            else:
                core.log_error("\n", "It is unclear how",
                               "\n", "the data stored in", key,
                               "\n", "should be concatenated...",
                               "\n", "For now, the data from all frames will",
                               "\n", "simply be copied into the combined frame.",
                               "\n", unit="concatdeconcat")
                self.copy_cartain_key_and_value_from_all_frames(
                    frames, key, combined_frame)
                combined_frame[key] = "ConcatenationMethodUnknown!"


        
        combined_frame["NeedsDeconcatenation"] = True
        return combined_frame
    
    
    def __call__(self, incoming_frame):
        
        if self.need_to_clear_cache:
            core.log_debug("\n", "Clearing the cache!"
                           "\n", unit="concatdeconcat")
            self.cache_of_frames        = [self.cache_of_frames[-1]]
            self.cache_of_frame_numbers = [self.cache_of_frame_numbers[-1]]
            self.need_to_clear_cache    = False
        
        if incoming_frame.type == core.G3FrameType.Scan:
            self.scan_frame_number += 1
            self.cache_of_frames.append(incoming_frame)
            self.cache_of_frame_numbers.append(self.scan_frame_number)
            core.log_debug("\n", "Added scan frame", self.scan_frame_number,
                           "to the cache.",
                           "\n", "Number of frames in cache:",
                           len(self.cache_of_frames),
                           "\n", "They are scan frames",
                           self.cache_of_frame_numbers,
                           "\n", unit="concatdeconcat")
            
            if len(self.cache_of_frames) < 3:
                return []
            elif ("Turnaround" not in self.cache_of_frames[0]) or \
                 ("Turnaround" not in self.cache_of_frames[2]) or \
                 ("Turnaround"     in self.cache_of_frames[1]):
                return self.cache_of_frames
            else:
                core.log_debug("\n", "The frames in the cache",
                               "are being combined...", "\n",
                               unit="concatdeconcat")
                combined_frame = self.concatenate_frames(self.cache_of_frames)
                core.log_debug("\n", "The contents of the combined frame:",
                               "\n", combined_frame,
                               "\n", unit="concatdeconcat")
                self.need_to_clear_cache = True
                return combined_frame
        
        if incoming_frame.type == core.G3FrameType.EndProcessing:
            if len(self.cache_of_frames) < 3:
                return self.cache_of_frames[1:] + [incoming_frame]
            else:
                combined_frame = self.concatenate_frames(self.cache_of_frames)
                self.need_to_clear_cache = True
                return combined_frame + [incoming_frame]

            


class DeconcatenateScanFrame(object):
    
    '''
    This pipeline module essentially does the inverse operations of
    what is done by the module ConcatenateScanFrames; it splits the data
    stored in a scan frame into three parts and returns three scan frames.
    
    '''
    
    def __init__(self, equivalent_keys={}):
        
        self.frame_counter   = 0
        self.equivalent_keys = equivalent_keys
    
    
    def deconcatenate_timestream(self, timestream, chunk_lengths):
        
        tss = []
        for i in range(len(chunk_lengths)):
            tss.append(timestream[sum(chunk_lengths[0:i]):
                                  sum(chunk_lengths[0:i+1])])
        return tss
    
    
    def deconcatenate_timestream_map(self, tsm, chunk_lens):
        
        tsms = [type(tsm)() for i in range(len(chunk_lens))]
        for k, ts in tsm.items():
            tss = self.deconcatenate_timestream(ts, chunk_lens)
            for i in range(len(chunk_lens)):
                tsms[i][k] = tss[i]
        return tsms
    
    
    def deconcatenate_g3_vector(self, vector, chunk_lengths):
        
        return self.deconcatenate_timestream(vector, chunk_lengths)
    
    
    def deconcatenate_map_g3_vector(self, mv, chnk_lens, n_frames=3):
        
        mvs = []
        for i in range(n_frames):
            mv_i = type(mv)()
            mvs.append(mv_i)
        for k, long_vctr in mv.items():
            chnks = chnk_lens[k]
            short_vctrs = self.deconcatenate_g3_vector(long_vctr, chnks)
            for i in range(len(short_vctrs)):
                mvs[i][k] = short_vctrs[i]
        return mvs
    
    
    def __call__(self, frame):
    
        if (frame.type == core.G3FrameType.Scan) and \
           ("NeedsDeconcatenation" in frame.keys()):
            
            core.log_debug("\n", "A scan frame that needs deconcatenation",
                           "has arrived.",
                           "\n", "The process is starting...",
                           "\n", unit="concatdeconcat")
            for key, value in self.equivalent_keys.items():
                core.log_debug("\n", key, "stores data that were newly added",
                               "\n", "by the processing module.",
                               "\n", "They will be deconcatenated in the",
                               "\n", "same way as the data stored in",
                               value+".", "\n", unit="concatdeconcat")
            
            self.frame_counter += 1
            
            list_of_frames = []
            n_frames = 3
            for i in range(n_frames):
                list_of_frames.append(core.G3Frame(core.G3FrameType.Scan))
            
            list_of_frames[0]["Turnaround"] = True
            list_of_frames[2]["Turnaround"] = True
            
            frame.pop("NeedsDeconcatenation")
            for key in frame.keys():
                if (key.startswith("(")):
                    pass
                
                elif key in ["GCPFeatureBits", "ObservationID", "SourceName"]:
                    for i in range(n_frames):
                        list_of_frames[i][key] = frame[key]
                
                elif key in ["DfMuxHousekeeping", "OnlinePointingModel",
                             "TrackerStatus"    , "TrackerPointing"]:
                    for i in range(n_frames):
                        suffix = "FromScanFrame"+str(i)+"InTheCache"
                        list_of_frames[i][key] = frame["("+key+suffix+")"]
                
                elif isinstance(frame[key], core.G3Timestream):
                    key_to_get_chnk_lens = \
                        "(" + self.equivalent_keys.get(key, key) + \
                        "ChunkLensUsedForConcatenation"+")"
                    chnk_lens = frame[key_to_get_chnk_lens]
                    tss = self.deconcatenate_timestream(frame[key], chnk_lens)
                    for i in range(n_frames):
                        list_of_frames[i][key] = tss[i]
                
                elif isinstance(frame[key], core.G3TimestreamMap):
                    key_to_get_chnk_lens = \
                        "(" + self.equivalent_keys.get(key, key) + \
                        "ChunkLensUsedForConcatenation"+")"
                    chnk_lens = frame[key_to_get_chnk_lens]
                    tsms = self.deconcatenate_timestream_map(
                               frame[key], chnk_lens)
                    for i in range(n_frames):
                        list_of_frames[i][key] = tsms[i]
                
                elif ("Vector"  in str(type(frame[key]))) and \
                     ("Map" not in str(type(frame[key]))):
                    key_to_get_chnk_lens = \
                        "(" + self.equivalent_keys.get(key, key) + \
                        "ChunkLensUsedForConcatenation"+")"
                    chnk_lens = frame[key_to_get_chnk_lens]
                    vctrs = self.deconcatenate_g3_vector(
                                frame[key], chnk_lens)
                    for i in range(n_frames):
                        list_of_frames[i][key] = vctrs[i]
                
                elif ("Vector"  in str(type(frame[key]))) and \
                     ("Map"     in str(type(frame[key]))):
                    key_to_get_chnk_lens = \
                        "(" + self.equivalent_keys.get(key, key) + \
                        "ChunkLensUsedForConcatenation"+")"
                    chnk_lens = frame[key_to_get_chnk_lens]
                    mvs = self.deconcatenate_map_g3_vector(
                              frame[key], chnk_lens)
                    for i in range(n_frames):
                        list_of_frames[i][key] = mvs[i]                
                
                else:
                    core.log_error("\n", "During the concatenation process,",
                                   "\n", "it was unclear what to do with",
                                   "\n", "the data stored in", key+",",
                                   "\n", "and the data from each frame were",
                                   "\n", "simply moved to the combined frame.",
                                   "\n", "Now they will be recovered.",
                                   "\n", unit="concatdeconcat")
                    for i in range(n_frames):
                        suffix = "FromScanFrame"+str(i)+"InTheCache"
                        list_of_frames[i][key] = frame["("+key+suffix+")"]                
            
            core.log_debug("\n", "Process complete!",
                           "\n", unit="concatdeconcat")
            
            if self.frame_counter == 1:
                return list_of_frames
            else:
                return list_of_frames[1:]




@core.pipesegment
def ConcatenateFramesAndProcessDataThenDeconcatenate(
        pipeline,
        processing_module=None,
        equivalent_keys={},
        **kwargs_for_processing):
    
    '''
    This pipeline segment concatenates three consecutive scan frames, which
    correspond to a turnaround, a constant speed scan, and another turnaround,
    applies certain user-specified processing step to the combined frame,
    and then deconcatenates the combined frame back into three scan frames.
    
    The original purpose of doing the processing on the combined frame is to
    reduce potential edge effects (such as oscillations) introduced by
    FFT-based timestream filters; since the data near the two ends of
    a timestream in this combined frame correspond to turnarounds,
    the potential edge effects can be mostly confined within these regions,
    and thus the data corresponding to the constant speed scan
    will be affected less.
    
    
    Arguments
    ---------
    processing_module: G3Pipeline module
        A pipeline module that will do some processing on the data
        stored in the combined frame.
    
    **kwargs_for_processing: keyword argements
        Keyword arguments such as input timestream name that are needed
        for the the processing module.
    
    equivalent_keys: dict
        If the processing module adds some output, say, a filtered
        timestream map, to the combined frame, this output needs to be
        properly deconcatenated. To specify how the output data should be
        deconcatenated (e.g. what length should be used for each of the
        three chunks into which timestreams are split),
        the user needs to provide a dictionary where, in a given
        key/value pair, the former is the key in the combined frame
        that stores the output data, and the latter is the key in the frame
        that stores some other data whose type is the same as the type of
        the output data. This way, whatever method is used to deconcatenate
        the latter data will be used to deconcatenate the output data.
        In the case of some timestream filter,
        {output_tsm_name: input_key_name} will probably work.
    
    '''
    
    pipeline.Add(ConcatenateScanFrames)
    
    if processing_module != None:
        pipeline.Add(processing_module, **kwargs_for_processing)

    pipeline.Add(DeconcatenateScanFrame,
                 equivalent_keys=equivalent_keys)





# Things related to time constant deconvolution
# -----------------------------------------------------------------------------

class TimeConstantDeconvolution(object):
    
    def __init__(self, input_tsm_name, output_tsm_name, tau_map_name):
        
        self.input_tsm_name  = input_tsm_name
        self.output_tsm_name = output_tsm_name
        self.tau_map_name    = tau_map_name
        self.tau_map         = None
    
    
    def __call__(self, frame):
        
        if self.tau_map_name in frame.keys():
            self.tau_map = frame[self.tau_map_name]
        
        if (frame.type == core.G3FrameType.Scan) and \
           (self.input_tsm_name in frame.keys()):
            if self.tau_map is None:
                frame[self.output_tsm_name] = frame[self.input_tsm_name]
            else:
                input_tsm    = frame[self.input_tsm_name]
                output_tsm   = core.G3TimestreamMap()
                dummy_filter = np.ones(input_tsm.n_samples)
                
                core.log_debug("\n", "Timestreams are being deconvolved",
                               "\n", "with the information on their taus...",
                               "\n", unit="taudeconvolution")
                todfilter.fft_filter_mem_friendly(
                    input_tsm, dummy_filter, output_tsm, True, self.tau_map)
                
                missing_bolos = set(input_tsm.keys()) - set(output_tsm.keys())
                for bolo in missing_bolos:
                    output_tsm[bolo] = input_tsm[bolo]
                
                frame[self.output_tsm_name] = output_tsm




def inject_fake_time_constants(frame):

    if "BolometerProperties" in frame.keys():
        tau_map = core.G3MapDouble()
        for bolo in frame["BolometerProperties"].keys():
            tau_map[bolo] = 5.0 * core.G3Units.ms
        frame["TimeConstants"] = tau_map





# Other things...
# -----------------------------------------------------------------------------

def add_new_key_to_scan_frames(frame):
    
    if frame.type == core.G3FrameType.Scan:
        frame["NewKey"] = "NewData"
        return



# =============================================================================





# =============================================================================
# Run a pipeline to test the modules above
# -----------------------------------------------------------------------------


# Define some global variables
# -----------------------------------------------------------------------------

core.set_log_level(core.G3LogLevel.LOG_DEBUG,
                   unit="concatdeconcat")

core.set_log_level(core.G3LogLevel.LOG_DEBUG,
                   unit="taudeconvolution")

random_files = ["/spt/data/bolodata/downsampled/ra0hdec-67.25/70907307/offline_calibration.g3",
                "/spt/data/bolodata/downsampled/ra0hdec-67.25/70907307/0000.g3"]



# Construct a pipeline and run it
# -----------------------------------------------------------------------------

pipeline = core.G3Pipeline()

pipeline.Add(core.G3Reader,
             filename=random_files,
             n_frames_to_read=12)

pipeline.Add(inject_fake_time_constants)

pipeline.Add(add_new_key_to_scan_frames)

pipeline.Add(ConcatenateFramesAndProcessDataThenDeconcatenate,
             processing_module=TimeConstantDeconvolution,
             input_tsm_name="RawTimestreams_I",
             output_tsm_name="DeconvolvedTimestreams_I",
             tau_map_name="TimeConstants",
             equivalent_keys={"DeconvolvedTimestreams_I": "RawTimestreams_I"})

pipeline.Add(core.Dump)

pipeline.Run(profile=True)


# =============================================================================

