#include <pybindings.h>
#include <string>

#include <mapmaker/aberrationmodules.h>
#include <mapmaker/aberrationutils.h>

void RelativisticAberrationPointing::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
    if (frame->type == G3Frame::Scan){
        auto alphas = frame->Get<G3MapVectorDouble>(detector_alpha_key_);
        auto deltas = frame->Get<G3MapVectorDouble>(detector_delta_key_);

		G3MapVectorDoublePtr alpha_out(new G3MapVectorDouble());
		G3MapVectorDoublePtr delta_out(new G3MapVectorDouble());

		aberrate_delta_alpha(*alphas, *deltas, coordinate_system_, *alpha_out,
                   *delta_out, peculiar_vel_, peculiar_ra_/G3Units::rad,
                   peculiar_dec_/G3Units::rad);
		frame->Put(detector_alpha_out_key_, alpha_out);
		frame->Put(detector_delta_out_key_, delta_out);
    }
    out.push_back(frame);
}

EXPORT_G3MODULE("mapmaker", RelativisticAberrationPointing,
        (init<double,double,double,MapCoordReference,std::string,std::string,
         std::string,std::string>(
            ( arg("peculiar_ra")="PeculiarRA",
              arg("peculiar_dec")="PeculiarDec",
              arg("peculiar_vel")="PeculiarVel",
              arg("coord_sys")="CoordinateSys",
			  arg("detector_alpha_key")="DetectorAlphaPointing",
			  arg("detector_delta_key")="DetectorDeltaPointing",
			  arg("detector_alpha_out_key")="OutputAlphaPointing",
			  arg("detector_delta_out_key")="OutputDeltaPointing"))),
        RELATIVISTIC_ABERRATION_POINTING_DOCSTR)
