#pragma once

#include <G3Module.h>
#include <G3Frame.h>
#include <coordinateutils/G3SkyMap.h>
#include <calibration/BoloProperties.h>

#define CALCULATE_DETECTOR_POINTING_DOCSTR                               \
        "Calculates the detectors offset from boresight in terms of ra/dec (or az/el) angle.\n\n" \
        "Alpha == Ra, Delta == Dec because we are physicists. \n" \
        "If timestream map key is set, uses the keys in that to decide which detectors to calculate things for.\n" \
        "If it is not present, uses the keys in bolo_properties to decide which  detectors to calculate things for.\n"  \
        "\n\n"                                   \
        "    bolo_properties (->BolometerPropertiesMap) : path to bolo geometry\n\n" \
        "    detector_transform_key (->G3VectorQuat): Transform that takes az=el=0 local coordinates to our boresight pointing in our desired coordinate system\n\n"                    \
        "    coord_sys (MapCoordRef): Fortunately for our telescope El and Dec are just a sign flip \n\n"                        \
        "    out_detector_alpha_key (->G3MapVectorDouble):  where to store the detector alpha pointing \n\n"                    \
        "    out_detector_delta_key (->G3MapVectorDouble):  where to store the detector delta pointing\n\n" \
        "    timestream_map_key (->G3TimestreamMap):  if not \"\", only calculate pointing for keys in the timestream map\n\n" \


class CalculateDetectorPointing : public G3Module 
{ 
public:
	CalculateDetectorPointing(
		std::string detector_transform_key,
		std::string bolo_props_key = "BolometerProperties",
		MapCoordReference coord_sys = Equatorial,
		std::string out_detector_alpha_key = "DetectorAlphaPointing",
		std::string out_detector_delta_key = "DetectorDeltaPointing",
		std::string timestream_map_key = "") :
		detector_transform_key_(detector_transform_key),
		bolo_props_key_(bolo_props_key),
		coord_sys_(coord_sys),
		out_detector_alpha_key_(out_detector_alpha_key),
		out_detector_delta_key_(out_detector_delta_key),
		timestream_map_key_(timestream_map_key)
	{}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:

	std::string detector_transform_key_;
	std::string bolo_props_key_;
	MapCoordReference coord_sys_;
	std::string out_detector_alpha_key_;
	std::string out_detector_delta_key_;
	std::string timestream_map_key_;
	BolometerPropertiesMapConstPtr bolo_props_;
};
