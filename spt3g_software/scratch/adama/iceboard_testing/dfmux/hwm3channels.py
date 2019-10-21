def info():
    channels = {'chan_1_1': 1013790.2,
    'chan_2_1': 1090822.0,
    'chan_3_1': 1201255.3,
    'chan_4_1': 1273024.0,
    'chan_5_1': 1387775.0,
    'chan_6_1': 1498169.5,
    'chan_7_1': 1584058.6,
    'chan_8_1': 1672748.4,
    'chan_9_1': 1738732.9,
    'chan_10_1': 1857610.3,
    'chan_11_1': 1958333.1,
    'chan_12_1': 2058316.5,
    'chan_13_1': 2157931.2,
    'chan_14_1': 2245638.4,
    'chan_15_1': 2324807.4,
    'chan_16_1': 2433113.1,
    'chan_17_1': 2545402.7,
    'chan_18_1': 2616332.3,
    'chan_19_1': 2701511.7,
    'chan_20_1': 2801930.1,
    'chan_21_1': 2903174.6,
    'chan_22_1': 2993206.7,
    'chan_23_1': 3104495.0,
    'chan_24_1': 3199664.8,
    'chan_25_1': 3301857.9,
    'chan_26_1': 3400870.9,
    'chan_27_1': 3482800.1,
    'chan_28_1': 3559158.7,
    'chan_29_1': 3680790.5,
    'chan_30_1': 3773860.8,
    'chan_31_1': 3880047.4,
    'chan_32_1': 3942431.6,
    'chan_33_1': 4055979.0,
    'chan_34_1': 4146607.1,
    'chan_35_1': 4220542.2,
    'chan_36_1': 4315948.0,
    'chan_37_1': 4416755.6,
    'chan_38_1': 4539549.4,
    'chan_39_1': 4597451.8,
    'chan_40_1': 4720183.1,
    'chan_41_1': 4823726.9,
    'chan_42_1': 4927233.5,
    'chan_43_1': 5001069.5,
    'chan_44_1': 5080165.2,
    'chan_45_1': 5178046.9,
    'chan_46_1': 5288286.4,
    'chan_47_1': 5359606.7,
    'chan_48_1': 5463401.1,
    'chan_49_1': 5575095.5,
    'chan_50_1': 5676031.1,
    'chan_51_1': 5771628.1,
    'chan_52_1': 5870601.3,
    'chan_53_1': 5945942.4,
    'chan_54_1': 6063321.8,
    'chan_55_1': 6161171.7,
    'chan_56_1': 6216051.6,
    'chan_57_1': 6352078.2,
    'chan_58_1': 6450803.5,
    'chan_59_1': 6514197.3,
    'chan_60_1': 6638750.7,
    'chan_61_1': 6713243.7,
    'chan_62_1': 6815352.0,
    'chan_63_1': 6912582.0,
    'chan_64_1': 6997413.8,
    'chan_1_2': 1013790.2,
    'chan_2_2': 1090822.0,
    'chan_3_2': 1201255.3,
    'chan_4_2': 1273024.0,
    'chan_5_2': 1387775.0,
    'chan_6_2': 1498169.5,
    'chan_7_2': 1584058.6,
    'chan_8_2': 1672748.4,
    'chan_9_2': 1738732.9,
    'chan_10_2': 1857610.3,
    'chan_11_2': 1958333.1,
    'chan_12_2': 2058316.5,
    'chan_13_2': 2157931.2,
    'chan_14_2': 2245638.4,
    'chan_15_2': 2324807.4,
    'chan_16_2': 2433113.1,
    'chan_17_2': 2545402.7,
    'chan_18_2': 2616332.3,
    'chan_19_2': 2701511.7,
    'chan_20_2': 2801930.1,
    'chan_21_2': 2903174.6,
    'chan_22_2': 2993206.7,
    'chan_23_2': 3104495.0,
    'chan_24_2': 3199664.8,
    'chan_25_2': 3301857.9,
    'chan_26_2': 3400870.9,
    'chan_27_2': 3482800.1,
    'chan_28_2': 3559158.7,
    'chan_29_2': 3680790.5,
    'chan_30_2': 3773860.8,
    'chan_31_2': 3880047.4,
    'chan_32_2': 3942431.6,
    'chan_33_2': 4055979.0,
    'chan_34_2': 4146607.1,
    'chan_35_2': 4220542.2,
    'chan_36_2': 4315948.0,
    'chan_37_2': 4416755.6,
    'chan_38_2': 4539549.4,
    'chan_39_2': 4597451.8,
    'chan_40_2': 4720183.1,
    'chan_41_2': 4823726.9,
    'chan_42_2': 4927233.5,
    'chan_43_2': 5001069.5,
    'chan_44_2': 5080165.2,
    'chan_45_2': 5178046.9,
    'chan_46_2': 5288286.4,
    'chan_47_2': 5359606.7,
    'chan_48_2': 5463401.1,
    'chan_49_2': 5575095.5,
    'chan_50_2': 5676031.1,
    'chan_51_2': 5771628.1,
    'chan_52_2': 5870601.3,
    'chan_53_2': 5945942.4,
    'chan_54_2': 6063321.8,
    'chan_55_2': 6161171.7,
    'chan_56_2': 6216051.6,
    'chan_57_2': 6352078.2,
    'chan_58_2': 6450803.5,
    'chan_59_2': 6514197.3,
    'chan_60_2': 6638750.7,
    'chan_61_2': 6713243.7,
    'chan_62_2': 6815352.0,
    'chan_63_2': 6912582.0,
    'chan_64_2': 6997413.8,
    'chan_1_3': 1013790.2,
    'chan_2_3': 1090822.0,
    'chan_3_3': 1201255.3,
    'chan_4_3': 1273024.0,
    'chan_5_3': 1387775.0,
    'chan_6_3': 1498169.5,
    'chan_7_3': 1584058.6,
    'chan_8_3': 1672748.4,
    'chan_9_3': 1738732.9,
    'chan_10_3': 1857610.3,
    'chan_11_3': 1958333.1,
    'chan_12_3': 2058316.5,
    'chan_13_3': 2157931.2,
    'chan_14_3': 2245638.4,
    'chan_15_3': 2324807.4,
    'chan_16_3': 2433113.1,
    'chan_17_3': 2545402.7,
    'chan_18_3': 2616332.3,
    'chan_19_3': 2701511.7,
    'chan_20_3': 2801930.1,
    'chan_21_3': 2903174.6,
    'chan_22_3': 2993206.7,
    'chan_23_3': 3104495.0,
    'chan_24_3': 3199664.8,
    'chan_25_3': 3301857.9,
    'chan_26_3': 3400870.9,
    'chan_27_3': 3482800.1,
    'chan_28_3': 3559158.7,
    'chan_29_3': 3680790.5,
    'chan_30_3': 3773860.8,
    'chan_31_3': 3880047.4,
    'chan_32_3': 3942431.6,
    'chan_33_3': 4055979.0,
    'chan_34_3': 4146607.1,
    'chan_35_3': 4220542.2,
    'chan_36_3': 4315948.0,
    'chan_37_3': 4416755.6,
    'chan_38_3': 4539549.4,
    'chan_39_3': 4597451.8,
    'chan_40_3': 4720183.1,
    'chan_41_3': 4823726.9,
    'chan_42_3': 4927233.5,
    'chan_43_3': 5001069.5,
    'chan_44_3': 5080165.2,
    'chan_45_3': 5178046.9,
    'chan_46_3': 5288286.4,
    'chan_47_3': 5359606.7,
    'chan_48_3': 5463401.1,
    'chan_49_3': 5575095.5,
    'chan_50_3': 5676031.1,
    'chan_51_3': 5771628.1,
    'chan_52_3': 5870601.3,
    'chan_53_3': 5945942.4,
    'chan_54_3': 6063321.8,
    'chan_55_3': 6161171.7,
    'chan_56_3': 6216051.6,
    'chan_57_3': 6352078.2,
    'chan_58_3': 6450803.5,
    'chan_59_3': 6514197.3,
    'chan_60_3': 6638750.7,
    'chan_61_3': 6713243.7,
    'chan_62_3': 6815352.0,
    'chan_63_3': 6912582.0,
    'chan_64_3': 6997413.8,
    'chan_1_4': 1013790.2,
    'chan_2_4': 1090822.0,
    'chan_3_4': 1201255.3,
    'chan_4_4': 1273024.0,
    'chan_5_4': 1387775.0,
    'chan_6_4': 1498169.5,
    'chan_7_4': 1584058.6,
    'chan_8_4': 1672748.4,
    'chan_9_4': 1738732.9,
    'chan_10_4': 1857610.3,
    'chan_11_4': 1958333.1,
    'chan_12_4': 2058316.5,
    'chan_13_4': 2157931.2,
    'chan_14_4': 2245638.4,
    'chan_15_4': 2324807.4,
    'chan_16_4': 2433113.1,
    'chan_17_4': 2545402.7,
    'chan_18_4': 2616332.3,
    'chan_19_4': 2701511.7,
    'chan_20_4': 2801930.1,
    'chan_21_4': 2903174.6,
    'chan_22_4': 2993206.7,
    'chan_23_4': 3104495.0,
    'chan_24_4': 3199664.8,
    'chan_25_4': 3301857.9,
    'chan_26_4': 3400870.9,
    'chan_27_4': 3482800.1,
    'chan_28_4': 3559158.7,
    'chan_29_4': 3680790.5,
    'chan_30_4': 3773860.8,
    'chan_31_4': 3880047.4,
    'chan_32_4': 3942431.6,
    'chan_33_4': 4055979.0,
    'chan_34_4': 4146607.1,
    'chan_35_4': 4220542.2,
    'chan_36_4': 4315948.0,
    'chan_37_4': 4416755.6,
    'chan_38_4': 4539549.4,
    'chan_39_4': 4597451.8,
    'chan_40_4': 4720183.1,
    'chan_41_4': 4823726.9,
    'chan_42_4': 4927233.5,
    'chan_43_4': 5001069.5,
    'chan_44_4': 5080165.2,
    'chan_45_4': 5178046.9,
    'chan_46_4': 5288286.4,
    'chan_47_4': 5359606.7,
    'chan_48_4': 5463401.1,
    'chan_49_4': 5575095.5,
    'chan_50_4': 5676031.1,
    'chan_51_4': 5771628.1,
    'chan_52_4': 5870601.3,
    'chan_53_4': 5945942.4,
    'chan_54_4': 6063321.8,
    'chan_55_4': 6161171.7,
    'chan_56_4': 6216051.6,
    'chan_57_4': 6352078.2,
    'chan_58_4': 6450803.5,
    'chan_59_4': 6514197.3,
    'chan_60_4': 6638750.7,
    'chan_61_4': 6713243.7,
    'chan_62_4': 6815352.0,
    'chan_63_4': 6912582.0,
    'chan_64_4': 6997413.8}
    return channels
