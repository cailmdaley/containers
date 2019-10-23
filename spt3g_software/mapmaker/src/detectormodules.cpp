#include <pybindings.h>
#include <mapmaker/detectormodules.h>
#include <coordinateutils/pointing.h>

void CalculateDetectorPointing::Process(G3FramePtr frame, std::deque<G3FramePtr> &out) 
{
     	if (frame->type == G3Frame::Calibration) {
		bolo_props_ = frame->Get<BolometerPropertiesMap>(bolo_props_key_);
	} else if (frame->type == G3Frame::Scan) { 
		if (!bolo_props_) {
			log_fatal("No bolometer properties seen before first data."
				  " Are you missing calibration data?");
		}
		G3VectorString dets;
		if (!timestream_map_key_.empty()){
			G3TimestreamMapConstPtr ts_map = 
			    frame->Get<G3TimestreamMap>(timestream_map_key_);
			dets = G3VectorString(ts_map->size());
			int i = 0;
			for (auto ts_id = ts_map->begin(); ts_id != ts_map->end(); ts_id++, i++){
				dets[i] = ts_id->first;
			}
		} else { 
			dets = G3VectorString(bolo_props_->size());
			int i = 0;
			for (auto ts_id=bolo_props_->begin(); ts_id!=bolo_props_->end(); 
			    ts_id++, i++){
				dets[i] = ts_id->first;
			}
		}

		G3VectorQuatConstPtr transforms = frame->Get<G3VectorQuat>(
		    detector_transform_key_);
		

		G3MapVectorDoublePtr alpha_out = boost::make_shared<G3MapVectorDouble>();
		G3MapVectorDoublePtr delta_out = boost::make_shared<G3MapVectorDouble>();
		for (auto it=dets.begin(); it!=dets.end(); it++){
			std::vector<double> alpha, delta;
			get_detector_pointing( bolo_props_->at(*it).x_offset,
					       bolo_props_->at(*it).y_offset,
					       *transforms, coord_sys_, 
					       alpha, delta);
			(*alpha_out)[*it] = alpha;
			(*delta_out)[*it] = delta;
		}
		frame->Put(out_detector_alpha_key_, alpha_out);
		frame->Put(out_detector_delta_key_, delta_out);
	}
	out.push_back(frame);		
}
namespace bp = boost::python;

EXPORT_G3MODULE("pointing", CalculateDetectorPointing,
		(init<std::string,std::string,MapCoordReference, std::string,
		 std::string,std::string>(
			(arg("detector_transform_key"),
			 arg("bolo_props_key") = "BolometerProperties",
			 arg("coord_sys")=Equatorial,
			 arg("out_detector_alpha_key") = "DetectorAlphaPointing",
			 arg("out_detector_delta_key") = "DetectorDeltaPointing",
			 arg("timestream_map_key") = ""))),
		CALCULATE_DETECTOR_POINTING_DOCSTR);
