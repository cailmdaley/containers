from spt3g import core
from spt3g.util import G3Viewer

f = core.G3Frame()
ts = core.G3Timestream([1,2,3,4])
f['A G3Timestream'] = ts
d = {'A Dict':
        {'Another Dict':
            {'A Third':
                {'A Final Dict':'WEEEE'}
            },
         'An Integer': 3,
        },
        'A List': [1, 2, 3, 4, [5, 6, 7, {'Wow, Another Dict!': 1}]],
        'A G3Frame': f
    }

v = G3Viewer(d)
