"""
This contains functionality to combine (add, subtract) frames of the same
Id in a G3 File and creates a new frame from that combined data
"""
from spt3g import core
from spt3g.util.genericutils import is_int
from spt3g.core import concatenate_timestreams
import copy
import fnmatch

def fc_add_vec_doub(a,b):
    return core.G3VectorDouble(np.asarray(a) + np.asarray(b))

def fc_add(a,b):
    return a + b

def fc_ident(a,b):
    return a

def fc_extend(a,b):
    a_cpy = copy.copy(a)
    a_cpy.extend(b)
    return a_cpy

def fc_add_int(a,b):
    return core.G3Int(core.G3Int(a).value + core.G3Int(b).value)

def fc_add_dic(a,b):
    '''
    Add the values in a dictionary
    '''
    b = copy.copy(b)
    for k in a.keys():
        if k in b:
            b[k] += a[k]
        else:
            b[k] = a[k]
    return b

def fc_concatenate_ts(a,b):
    return concatenate_timestreams([a,b])

@core.indexmod
class FrameCombiner(object):
    '''
    Frame Combiner is a class for combining data from multiple frames into 
    one single frame. It is only meant to work on one type of frame specified 
    by 'type'
    
    Right now it collects all of the frames until the end of processing.
    If fr_id is set, it will only collect frames with a matching 'Id' key.

    It only collects the keys specified in key_op_mapping.
    It performs some form of addition on them.

    It outputs one frame with all the of the specified keys "added".
    If fr_id is set, the output frame will have the id `'combined_'+fr_id`

    key_op_mapping is a dictionary of the form `{ Frame Key: operation }`

    operation is a function that accepts two arguments.  
    It returns some form of combination/summing
    on the two arguments, adding, appending, etc.  

    
    '''

    def __init__(self, type, key_op_mapping, fr_id=None):
        self.frame_type = type
        self.output_frame = core.G3Frame(type)
        self.kop = key_op_mapping
        self.id = fr_id
        if self.id:
            self.output_frame['Id'] = 'combined_'+self.id
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            if list(self.output_frame.keys()) in [ [], ['Id'] ]:
                return
            else:
                return [self.output_frame, frame]
        if frame.type != self.frame_type:
            return
        if self.id:
            if not ('Id' in frame and fnmatch.fnmatch(frame['Id'], self.id)):
                return
        for k, op in self.kop.items():
            if not k in frame:
                continue
            if not k in self.output_frame:
                try:
                    self.output_frame[k] = frame[k]
                except:
                    if is_int(frame[k]):
                        self.output_frame[k] = core.G3Int(frame[k])
            else:
                oframe_val = self.output_frame[k]
                del self.output_frame[k]
                try:
                    self.output_frame[k] = op(oframe_val, frame[k])
                except:
                    print(type(oframe_val), oframe_val)
                    print(k)
                    raise RuntimeError()
        return []

@core.indexmod
class MapFrameCombiner(object):
    def __init__(self, fr_id=None):
        self.fc = FrameCombiner(type = core.G3FrameType.Map,
                                fr_id = fr_id,
                                key_op_mapping = {
                'T' : fc_add,
                'Q' : fc_add,
                'U' : fc_add,
                'Wpol' : fc_add,
                'Wunpol' : fc_add})
    def __call__(self, frame):
        return self.fc(frame)
