from spt3g import core, dfmux, calibration
import numpy

# Utilities to invert crosstalk from detector timestreams

def invert_xtalk_matrix(xtalkcomponents):
    '''
    Invert a sparse crosstalk matrix. Takes (and returns) a double dictionary
    of the form xtalk['receiving_bolo']['source_bolo'] = fraction_of_source_bolo
    '''
    bolos = xtalkcomponents.keys()
    xtalkmat = numpy.matrix(numpy.eye(len(bolos)))

    # Form dense matrix from sparse storage
    for i in range(len(bolos)):
        bolo = bolos[i]
        
        for component in xtalkcomponents[bolo].keys():
            if component not in bolos:
                continue
            val = xtalkcomponents[bolo][component]
            xtalkmat[i, bolos.index(component)] = val

    # Invert matrix
    xtalkmat_inv = numpy.linalg.inv(xtalkmat)

    # Sparsify and index output
    xtalkcomponents_inv = {}

    for i in range(len(bolos)):
        bolo = bolos[i]
        xtalkcomponents_inv[bolo] = {}
        for component in xtalkcomponents[bolo].keys():
            xtalkcomponents_inv[bolo][component] = \
              float(xtalkmat_inv[i, bolos.index(component)])
            # Conversion to float is to avoid slowdowns with 0-dim arrays
    return xtalkcomponents_inv

@core.indexmod
class InvertCrosstalk(object):
    '''
    Invert crosstalk from timestreams by applying the inverse of the crosstalk
    matrix stored in the calibration frame under the key "CrosstalkMatrix".
    Input timestreams can be in any units, though the efficiency of the
    calculation will be highest if you provide them in Current first (otherwise
    they are internally converted to Current). Output is always in Current
    units if a calculation is done (see note about ignore_missing).
    If 'ignore_missing' is set to False (the default), will
    assert if the crosstalk matrix is not present; otherwise, will just copy
    its input to its output, as though the matrix were the identity matrix.
    '''

    def __init__(self, input='Timestreams',
      output='CrosstalkCleanedTimestreamsAmps', ignore_missing=False):
        self.input = input
        self.output = output
        self.matrix = None
        self.ignore_missing = ignore_missing
        self.captive_unit_converter_out = 'XXXTEMPAMPSXXX_' + input
        self.captive_unit_converter = dfmux.ConvertTimestreamUnits(Input=input,
            Units=core.G3TimestreamUnits.Current,
            Output=self.captive_unit_converter_out)

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration and \
          'CrosstalkMatrix' in frame:
            self.matrix = invert_xtalk_matrix(frame['CrosstalkMatrix'])

        if frame.type != core.G3FrameType.Scan:
           self.captive_unit_converter(frame)
           return 

        input = frame[self.input]

        # Do null inversion if asked
        if self.matrix is None and self.ignore_missing:
            frame[self.output] = input
            return

        if input.values()[0].units != core.G3TimestreamUnits.Current:
            self.captive_unit_converter(frame)
            input = frame.pop(self.captive_unit_converter_out)

        # Record which timestreams have bad data and should not be mixed
        badts = []
        for xbolo, xts in input.items():
            if not numpy.isfinite(xts).all():
                badts.append(xbolo)

        output = core.G3TimestreamMap()
        for bolo,ts in input.iteritems():
            if bolo not in self.matrix.keys() or bolo in badts:
                output[bolo] = ts
                continue

            newts = core.G3Timestream(ts)
            newts[:] = numpy.zeros(len(newts))
            
            # Output TS = matrix * timestreams
            for xbolo,xtalk in self.matrix[bolo].items():
                if xbolo not in input or xbolo in badts:
                    continue
                xts = input[xbolo]
                newts._assert_congruence(xts)
                newts += xtalk*numpy.asarray(xts)
                # numpy.asarray here avoids a G3Timestream copy constructor,
                # which speeds up processing. It also skips some error
                # checking, which is less ideal, so that is done above.

            output[bolo] = newts

        frame[self.output] = output

@core.pipesegment
def CrosstalkCleanedTimestreams(pipe, input='RawTimestreams_I',
                                output='CrosstalkCleanedTimestreams_I', 
                                units=core.G3TimestreamUnits.Power,
                                ignore_missing=False):
    '''
    Create crosstalk-cleaned timestreams in the desired units from a set of
    input timestreams in arbitray other units. Handles the appropriate unit
    conversions for crosstalk inversion.
    '''
    pipe.Add(InvertCrosstalk, input=input, ignore_missing=ignore_missing,
             output='XXXTEMPAMPSCLEANXXX_' + output)
    pipe.Add(dfmux.ConvertTimestreamUnits, Input='XXXTEMPAMPSCLEANXXX_' + output, 
             Output=output, Units=units)
    pipe.Add(core.Delete, keys=['XXXTEMPAMPSCLEANXXX_' + output])

@core.indexmod
def LoadTextCrosstalkMatrix(frame, matrixpath=None, ignore_missing = False):
    '''
    Load a legacy (SPTpol-style) text crosstalk matrix into the calibration
    frame
    '''

    if frame.type != core.G3FrameType.Calibration:
        return

    xtalk = core.G3MapMapDouble()

    textmat = open(matrixpath, 'r')
    lines = textmat.readlines()
    textmat.close()

    # Form mapping from physical ID (used for SPTpol matrices) to logical ID
    bpm = frame['BolometerProperties']
    id_mapping = {}
    for bolo in bpm.keys():
        id_mapping[bpm[bolo].physical_name] = bolo

    # First line is supposed to be number of lines
    assert(len(lines) == int(lines[0]) + 1)
    lines = lines[1:]

    # Parse file
    for line in lines:
        fields = line.split() # RX_bolo TX_bolo xtalk_fraction

        if ignore_missing:
            if ((not fields[0] in id_mapping) or 
                (not fields[1] in id_mapping)):
                continue
        
        # Map physical to logical IDs
        fields[0] = id_mapping[fields[0]]
        fields[1] = id_mapping[fields[1]]

        if fields[0] not in xtalk:
            xtalk[fields[0]] = core.G3MapDouble()

        # Force crosstalk diagonal to be equal to 1 as SPTpol did
        if fields[0] != fields[1]:
            xtalk[fields[0]][fields[1]] = float(fields[2])
        else:
            xtalk[fields[0]][fields[1]] = 1.

    frame['CrosstalkMatrix'] = xtalk

