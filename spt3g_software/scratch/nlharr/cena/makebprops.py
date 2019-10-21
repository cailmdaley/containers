from spt3g import core, calibration
cena_line_up = {
    'cena-20150312_131911.g3':('rcw38-20150312_123711.g3', 'calib-20150311_135102.g3'),
    'cena-20150325_151459.g3':('rcw38-20150325_143300.g3', 'calib-20150325_143143.g3'),
    'cena-20150301_230758.g3':('rcw38-20150301_222600.g3', 'calib-20150301_222443.g3'),
    'cena-20150314_025927.g3':('rcw38-20150314_021727.g3', 'calib-20150314_021611.g3'),
    'cena-20150303_053730.g3':('rcw38-20150303_045532.g3', 'calib-20150303_045415.g3'),
    'cena-20150315_161543.g3':('rcw38-20150315_153345.g3', 'calib-20150315_153229.g3'),
    'cena-20150303_122107.g3':('rcw38-20150303_113907.g3', 'calib-20150303_113751.g3'),
    'cena-20150317_065622.g3':('rcw38-20150317_061426.g3', 'calib-20150317_061309.g3'),
    'cena-20150304_174541.g3':('rcw38-20150304_170344.g3', 'calib-20150304_170227.g3'),
    'cena-20150318_194938.g3':('rcw38-20150318_190454.g3', 'calib-20150318_190337.g3'),
    'cena-20150401_014312.g3':('rcw38-20150401_103312.g3', 'calib-20150401_100411.g3'),
    'cena-20150306_071532.g3':('rcw38-20150306_063336.g3', 'calib-20150307_063009.g3'),
    'cena-20150320_091251.g3':('rcw38-20150320_083055.g3', 'calib-20150320_082938.g3'),
    'cena-20150307_204739.g3':('rcw38-20150307_200253.g3', 'calib-20150307_200137.g3'),
    'cena-20150321_225618.g3':('rcw38-20150321_221418.g3', 'calib-20150321_221302.g3'),
    'cena-20150309_094036.g3':('rcw38-20150309_085837.g3', 'calib-20150309_085721.g3'),
    'cena-20150322_114432.g3':('rcw38-20150322_110232.g3', 'calib-20150322_114316.g3'),
    'cena-20150310_235454.g3':('rcw38-20150310_231255.g3', 'calib-20150310_231139.g3'),
    'cena-20150324_013009.g3':('rcw38-20150324_004811.g3', 'calib-20150324_004655.g3')
}

calib_dir = '/data52/nwhitehorn/3g-calib/calibrator/results/'
rcw38_cal_dir = '/data52/nwhitehorn/3g-calib/rcw38/results/'
out_dir = '/data/nlharr/cena_maps/cal_frames/'

for k, v in cena_line_up.items():
    in_fn = rcw38_cal_dir+v[0]
    out_fn = out_dir + k

    pipe = core.G3Pipeline()

    pipe.Add(core.G3Reader, filename=core.G3VectorString([rcw38_cal_dir+v[0], calib_dir + v[1]]))
    pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration or fr.type == core.G3FrameType.EndProcessing)


    pipe.Add(calibration.build_cal_frames.BuildBoloPropertiesMap, 
             drop_original_frames=False, filter_abs_point=True)
    pipe.Add(calibration.build_cal_frames.MergeCalibrationFrames, KeysToIgnore=[])

    pipe.Add(lambda fr: fr.type == core.G3FrameType.Calibration)
    pipe.Add(core.Dump)


    pipe.Add(core.G3Writer, filename=out_fn)
    pipe.Run()
