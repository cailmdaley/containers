#pragma once
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Map.h>
#include <G3Vector.h>
#include <G3Timestream.h>
#include <G3Data.h>
#include <G3Logging.h>

#include <coordinateutils/G3SkyMap.h>
#include <calibration/BoloProperties.h>


#define MAP_POINTING_CALCULATOR_DOCSTR                                  \
        "Figures out which bin each bolometer sample lives in.\n\n"     \
        "Does this for every key in the G3TimestreamMap at ts_map_key"  \
        "\n\n"                                                          \
        "    map_id [std::string] The id used to identify the map to calculate pointing for \n\n" \
        "    ts_map_key [->G3TimestreamMap] pointing is only computed for the keys in ts_map \n\n"  \
        "    trans_key (->G3VectorQuat): Transform that takes az=el=0 local coordinates to our boresight pointing in our desired coordinate system\n\n" \
        "    bolo_props_key [->BolometerPropertiesMap]  \n\n"  \
        "    detector_alpha_key [->G3MapVectorDouble]: If we've precomputed the" \
	" pointing this is where the alpha angle lives.  Usually you want this to be \"\" \n\n" \
        "    detector_delta_key [->G3MapVectorDouble]: If we've precomputed the" \
	" pointing this is where the delta angle lives.  Usually you want this to be \"\" \n\n" \
        "    pointing_store_location [->G3MapVectorInt]:  The map bin the bolo timestream sample is stored in.\n\n"

class MapPointingCalculator : public G3Module{
public:
	MapPointingCalculator(std::string map_id, 
			       std::string ts_map_key,

			       std::string bolo_props_key,
			       std::string trans_key,
			       
			       std::string detector_alpha_key,
			       std::string detector_delta_key,
			       
			       std::string pointing_store_key) :
		map_id_(map_id),ts_map_key_(ts_map_key),
		store_loc_(pointing_store_key), 
		bolo_props_key_(bolo_props_key),
		trans_key_(trans_key),
		detector_alpha_key_(detector_alpha_key), 
		detector_delta_key_(detector_delta_key) {
		pre_computed_pointing_ = detector_alpha_key != "" && detector_delta_key != "";
	}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string map_id_;
	std::string ts_map_key_;
	std::string store_loc_;

	std::string bolo_props_key_;
	std::string trans_key_;
	std::string detector_alpha_key_;
	std::string detector_delta_key_;

	boost::shared_ptr< const G3SkyMap > map_;
	BolometerPropertiesMapConstPtr bolo_props_;

	bool pre_computed_pointing_;

	SET_LOGGER("MapPointingCalculator");
};

#define MAP_BINNER_DOCSTR   \
        "Adds the timestream data to the maps present in the frame with 'Id' == map_id"\
        " It will try to add weight/timestream data for all of the maps present in the map frame specified by Id"\
        " By that I mean, it will only add weight information, or Q/U information if those maps are present in the input map frame.\n\n"\
        "    map_id [std::string]: id of map frame to use\n\n"  \
        "    ts_map_key [->G3TimestreamMap]: Timestreams you want to bin\n\n"\
        "    pointing_key [->G3MapVectorInt]: Where each timestream is pointing in the map\n\n"\
        "    weight_key [->G3MapDouble]: The weight assigned to each timestream\n\n"    \
        "    bolo_props_key [->BolometerPropertiesMap]  \n\n"  \
        "    bolo_ids_to_map_key [->G3String]: If \"\", use the keys in the ts_map.  If it is not \"\", maps to a list of detectors to actually bin\n\n"\
        "    trans_key (->G3VectorQuat): Transform that takes az=el=0 local" \
	"coordinates to our boresight pointing in our desired coordinate " \
	"system.  You only need to specify this if you are including the " \
	"polarization rotation.\n\n" \
        "    include_pol_rotation (bool):  If true accounts for the rotation of the detector polarization from the coordinate transform\n\n" 


class MapBinner : public G3Module{
public:
	MapBinner(std::string map_id, 
		   std::string ts_map_key, 
		   std::string pointing_key,  
		   std::string weight_key, 
		   std::vector<std::string> bolo_ids_to_map,
		   std::string bolo_props_key,
		   std::string trans_key,
		   bool include_pol_rotation
		) :
		map_id_(map_id), weight_key_(weight_key), 
		ts_map_key_(ts_map_key), pointing_key_(pointing_key),
		bolo_props_key_(bolo_props_key),
		bolo_ids_to_map_(bolo_ids_to_map), trans_key_(trans_key),
		include_pol_rotation_(include_pol_rotation)
	{}	
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string map_id_;
	bool is_polarized_;
	G3FramePtr map_frame_;

	std::string weight_key_;
	std::string ts_map_key_;
	std::string pointing_key_;
	std::string bolo_props_key_;

	std::vector<std::string> bolo_ids_to_map_;

	bool is_pol_;
	bool is_weight_;

	G3SkyMapPtr T,Q,U;
	G3SkyMapWeightsPtr W;
	BolometerPropertiesMapConstPtr bolo_props_;

	std::string trans_key_;
	bool include_pol_rotation_;

	SET_LOGGER("MapBinner");
};



class SimulatedTimestreamFiller : public G3Module {
public:
	SimulatedTimestreamFiller(std::string map_id,
				  std::string sim_pointing_key, 
				  std::string bolo_props_key, 
				  std::string ts_to_get_sample_rate_key,
				  std::string trans_key,
				  std::string ts_lst_key,
				  bool include_pol_rot,
				  std::string out_ts_key) :
		map_id_(map_id), pointing_key_(sim_pointing_key), 
		bolo_props_key_(bolo_props_key),
		ts_to_get_sample_rate_key_(ts_to_get_sample_rate_key),
		trans_key_(trans_key),
		include_pol_rot_(include_pol_rot),
		out_ts_key_(out_ts_key),
		ts_lst_key_(ts_lst_key)
	{}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string map_id_;
	std::string pointing_key_;
	std::string bolo_props_key_;
	std::string ts_to_get_sample_rate_key_;
	std::string trans_key_;
	bool include_pol_rot_;
	std::string out_ts_key_;
	std::string ts_lst_key_;

	BolometerPropertiesMapConstPtr bolo_props_;
	G3SkyMapConstPtr T,Q,U;
	bool is_pol_;
};




class SkyMapInterpTodFiller : public G3Module {
public:
	SkyMapInterpTodFiller(std::string map_id,
				std::string bolo_props_key, 
				std::string ts_to_get_sample_rate_key,
				std::string trans_key,
				bool include_pol_rot,
				std::string valid_ids,
				std::string detector_alpha_key,
				std::string detector_delta_key,
				std::string out_ts_key) :
		map_id_(map_id), 
		bolo_props_key_(bolo_props_key),
		ts_to_get_sample_rate_key_(ts_to_get_sample_rate_key),
		trans_key_(trans_key),
		include_pol_rot_(include_pol_rot),
		valid_ids_(valid_ids),
		detector_alpha_key_(detector_alpha_key),
		detector_delta_key_(detector_delta_key),
		out_ts_key_(out_ts_key)
	{pre_computed_pointing_ = detector_alpha_key != "" && detector_delta_key != "";}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string map_id_;
	std::string bolo_props_key_;
	std::string ts_to_get_sample_rate_key_;
	std::string trans_key_;
	bool include_pol_rot_;
	std::string valid_ids_;
	std::string detector_alpha_key_;
	std::string detector_delta_key_;
	std::string out_ts_key_;
	BolometerPropertiesMapConstPtr bolo_props_;
	G3SkyMapConstPtr T,Q,U;
	bool pre_computed_pointing_;
	SET_LOGGER("SimulatedTimestreamFiller");
};
