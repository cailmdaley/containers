from spt3g import core

@core.indexmod
class AddMetaData(object):
    '''
    This class will record the map-making parameters used and collect
    statistics on live- and dropped-bolometers in each scan. This information
    is stored in a PipelineInfo frame output with the rest of the data.

    Input variables:

    stats: an instantiation of the "timestreamflagging.GenerateFlagStats" class,
      to have been inserted into the pipeline before AddMetaData is called.

    args: a dictionary containing the script's input arguments.  The key should
      tell what the argument is (i.e. 'pixel_resolution'), and the value is,
      clearly, the argument's value.

    Example usage:
      parser = argparse.ArgumentParser()
      <add arguments>
      args = parser.parse_args()
      <do pipeline things, including flagging bolometers>
      stats = timestreamflagging.GenerateFlagStats(flag_key='Flags')
      pipe.Add(stats)
      <other pipeline modules>
      pipe.Add(std_processing.AddMetaData, stats = stats, args = vars(args))

      Note the use of vars(args) to turn the argparse object into a dictionary.
    '''
    def __init__(self, stats, args):
        self.bolo_list = core.G3MapVectorString()
        self.scan_number = 0
        self.args = args
        self.stats = stats

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan:
            # Record bolos surviving flagging
            if 'DeflaggedTimestreams' in frame:
                self.bolo_list[str(self.scan_number).zfill(3)] = frame[
                  'DeflaggedTimestreams'].keys()
            self.scan_number += 1

        if frame.type == core.G3FrameType.EndProcessing:
            meta = core.G3Frame(core.G3FrameType.PipelineInfo)
            # Record all mapmaking parameters
            for k, v in self.args.items():
                if v is None:
                    v = ''
                if isinstance(v, list):
                    v = core.G3VectorString([str(x) for x in v])
                meta[k] = v
            meta['SurvivingBolos'] = self.bolo_list
            # Record which bolos received which flags
            flag_stats = core.G3MapVectorVectorString()
            for scan_ind, flagd in enumerate(self.stats.scan_bolos):
                for reason, bolos in flagd.items():
                    if reason not in flag_stats:
                        flag_stats[reason] = [core.G3VectorString() for
                          x in range(scan_ind)]
                    if len(flag_stats[reason]) < scan_ind:
                        flag_stats[reason].extend([core.G3VectorString()
                          for x in range(scan_ind-len(flag_stats[reason]))])
                    flag_stats[reason].append(bolos)
            for reason, arr in flag_stats.items():
                if len(arr)<scan_ind+1:
                    flag_stats[reason].extend([core.G3VectorString()
                      for x in range(scan_ind+1-len(arr))])
            meta['DroppedBolos'] = flag_stats

            return [meta, frame]
