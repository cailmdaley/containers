#pragma once
#include <string>
#include <deque>
#include <G3Frame.h>
#include <G3Module.h>
#include <G3Logging.h>

#include <coordinateutils/G3SkyMap.h>

class RelativisticAberrationPointing : public G3Module{
public:
	RelativisticAberrationPointing(double peculiar_ra,
	                               double peculiar_dec,
	                               double peculiar_vel,
	                               MapCoordReference coord_sys,
	                               std::string detector_alpha_key,
	                               std::string detector_delta_key,
	                               std::string detector_alpha_out_key,
	                               std::string detector_delta_out_key) :
		peculiar_ra_(peculiar_ra), peculiar_dec_(peculiar_dec),
		peculiar_vel_(peculiar_vel), coordinate_system_(coord_sys),
		detector_alpha_key_(detector_alpha_key),
		detector_delta_key_(detector_delta_key),
		detector_alpha_out_key_(detector_alpha_out_key),
		detector_delta_out_key_(detector_delta_out_key){}
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	double peculiar_ra_;
	double peculiar_dec_;
	double peculiar_vel_;
	MapCoordReference coordinate_system_;
	std::string detector_alpha_key_;
	std::string detector_delta_key_;
	std::string detector_alpha_out_key_;
	std::string detector_delta_out_key_;
	SET_LOGGER("RelativisticAberrationPointing");
};

#define RELATIVISTIC_ABERRATION_POINTING_DOCSTR   \
	"Adds pointing offsets due to relativistic aberration\n\n" \
	"   peculiar_ra [double]: Right ascension of the observer's peculiar velocity\n\n" \
	"   peculiar_dec [double]: Declination of the observer's peculiar velocity\n\n" \
	"   peculiar_vel [double]: Magnitude of observer's peculiar velocity\n\n" \
	"   detector_alpha_key [->G3MapVectorDouble]:  Where the detector alpha pointings is stored\n\n" \
	"   detector_delta_key [->G3MapVectorDouble]: Where the detector delta pointing is stored\n\n" \
	"   detector_alpha_out_key [->G3MapVectorDouble]:  Where the aberrated alpha pointings is to be stored\n\n" \
	"   detector_delta_out_key [->G3MapVectorDouble]: Where the aberrated delta pointing is to be stored\n\n"
