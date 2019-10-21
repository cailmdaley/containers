import argparse
import pickle
from spt3g import core, gcp
import nhutils as nhu

def get_T(arcfiles, outfile):
    prefix = ['array', 'cryo', 'temperature', 0]
    keys = {'uchead': prefix + [0],
            'ichead': prefix + [1],  
            'ucstage': prefix + [10],
            'lctower': prefix + [11],
            'icstage': prefix + [12],
            '4khead': prefix + [13], 
            '50khead': prefix + [15],
            'time': ['array', 'cryo', 'utc'],
            'az': ['antenna0', 'tracker', 'actual', 0],
            'el': ['antenna0', 'tracker', 'actual', 1],
            'optics_4khead': ['array', 'cryo', 'temperature', 2, 8],
            'optics_50khead': ['array', 'cryo', 'temperature', 2, 4],
            }
    
    out = nhu.extract_keys(arcfiles, keys, frametype = core.G3FrameType.GcpSlow)
    f = open(outfile, 'wb')
    pickle.dump(out, f)
    f.close()
    return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', nargs = '+', type = str)
    parser.add_argument('-o', '--output', type = str, default = 'out.pkl')
    args = parser.parse_args()
    get_T(args.input, args.output)
