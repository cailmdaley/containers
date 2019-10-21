#include <pybindings.h>
#include <G3.h>
#include <G3Module.h>
#include <G3Logging.h>

#include <string>
#include <mapmaker/mappingmodules.h>
#include <coordinateutils/pointing.h>
#include <mapmaker/todmapconverters.h>

void MapPointingCalculator::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if (frame->type == G3Frame::Map){
		log_debug("running map");
		if (frame->Get<G3String>("Id")->value == map_id_){
			g3_assert(frame->Has("T"));
			map_ = frame->Get<G3SkyMap>("T");
		}
	} else if (frame->type == G3Frame::Calibration && frame->Has(bolo_props_key_)){
		bolo_props_ = frame->Get<BolometerPropertiesMap>(bolo_props_key_);
	} else if (frame->type == G3Frame::Scan) {
		if (! map_) {
			log_fatal("MapPointingCalculator did not receive a map frame with"
				  "specified id\n");
		}

		G3MapVectorIntPtr output_indices = 
			boost::make_shared<G3MapVectorInt>();
		G3TimestreamMapConstPtr tsm = 
			frame->Get<G3TimestreamMap>(ts_map_key_);
		std::vector<std::string> dets;
		for (auto it=tsm->begin(); it!=tsm->end(); it++){
			dets.push_back(it->first);
			(*output_indices)[it->first] = 
				std::vector<int>(it->second->size());
		}
		if (pre_computed_pointing_){
			G3MapVectorDoubleConstPtr alpha_scan_det =
				frame->Get<G3MapVectorDouble>(detector_alpha_key_);
			G3MapVectorDoubleConstPtr delta_scan_det = 
				frame->Get<G3MapVectorDouble>(detector_delta_key_);
			size_t i;
			for (i=0; i < dets.size(); i++){
				(*output_indices)[dets[i]] = map_->angles_to_pixels(
					alpha_scan_det->at(dets[i]),
					delta_scan_det->at(dets[i]));
			}
			frame->Put(store_loc_, output_indices);
		} else {
			G3VectorQuatConstPtr trans = frame->Get<G3VectorQuat>(trans_key_);
			size_t i;

			std::vector<double> alpha(trans->size());
			std::vector<double> delta(trans->size());

			for (i=0; i < dets.size(); i++){
				if (bolo_props_->count(dets[i])){
					const BolometerProperties & bp = bolo_props_->at(dets[i]);
					get_detector_pointing(bp.x_offset, bp.y_offset,
							      *trans, map_->coord_ref, alpha, delta);
					std::vector<int> & ref = (*output_indices)[dets[i]];
					ref = map_->angles_to_pixels(alpha, delta);
				} else{
					log_fatal("Bolo %s missing entry in BolometerPropertiesMap; "
						  "cannot calculate pointing.\n", dets[i].c_str());
				}
			}

			frame->Put(store_loc_, output_indices);
		}
	}
	out.push_back(frame);
}

void MapBinner::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	log_debug("Starting Process");
	if (frame->type == G3Frame::Scan){
		// Ensure a map is allocated
		if (!map_frame_) {
			log_fatal("No map frame");
		}

		// Skip frames that may not have data from our detectors in this
		// scan
		if (!frame->Has(ts_map_key_)){
			log_warn("Timestream Key %s is missing from the Scan "
				 "frame data for map binning.", ts_map_key_.c_str());
			out.push_back(frame);
			return;
		}

		// Trivial binning for a single timestream
		if (!is_pol_ && frame->Has<G3Timestream>(ts_map_key_)) {
			G3TimestreamConstPtr ts = frame->Get<G3Timestream>(ts_map_key_);
			G3VectorIntConstPtr pnt = frame->Get<G3VectorInt>(pointing_key_);
			g3_assert(pnt);
			if (T->units == G3Timestream::None)
				T->units = ts->units;
			if (T->units != ts->units)
				log_fatal("attempting to make maps with multiple unit types");
			bin_tod_pol(*ts,0,0,*pnt,1.0,*T,*T,*T);
			if (is_weight_){
				bin_w(*ts,*pnt,1.0,*W);
			}
			out.push_back(frame);
			return;
		}

		// Validate frame data
		if (!bolo_props_){
			log_fatal("no bolo props");
		}
		g3_assert(frame->Has<G3TimestreamMap>(ts_map_key_));		
		g3_assert(frame->Has<G3MapDouble>(weight_key_));
		g3_assert(frame->Has<G3MapVectorInt>(pointing_key_));

		if (frame->Get<G3TimestreamMap>(ts_map_key_)->begin()->second->size() !=
		    frame->Get<G3MapVectorInt>(pointing_key_)->begin()->second.size() ) {
			log_fatal("Pointing and detector timestream lengths differ.");
		}
		//grab the info we want
		G3MapDoubleConstPtr detector_weights = frame->Get<G3MapDouble>(weight_key_);
		boost::shared_ptr<const G3MapVectorInt> det_pnting = 
			frame->Get<G3MapVectorInt>(pointing_key_);
		G3TimestreamMapConstPtr timestreams = frame->Get<G3TimestreamMap>(ts_map_key_);
		//figure out which detectors we want info from

		std::vector<std::string> ts_ids;
		if (! bolo_ids_to_map_.empty()){
			ts_ids = bolo_ids_to_map_;
		} else{
			for (auto iter = timestreams->begin(); 
			     iter != timestreams->end(); iter++){
				ts_ids.push_back(iter->first);
			}
		}

		if (T->units == G3Timestream::None && timestreams->size() > 0){
			G3Timestream::TimestreamUnits u = timestreams->begin()->second->units;
			T->units = u;
			if (is_pol_){
				Q->units = u;
				U->units = u;
			}
		} 

		if (T->units != G3Timestream::None && 
		    T->units != timestreams->begin()->second->units){
			log_fatal("attempting to make maps with multiple unit types");
		}


		//actually bin the TOD
		//unpolarized caes
		if (! is_pol_){
			for (auto it=ts_ids.begin(); it!=ts_ids.end();it++){
				if (timestreams->find(*it) == timestreams->end()){
					log_info("Skipping %s", it->c_str());
					continue;					
				}
				if (det_pnting->find(*it) == det_pnting->end()){
                                        log_fatal("Unable to find pointing information for %s",
                                                  it->c_str());
                                }
				const std::vector<int>pnt=det_pnting->at(*it);
				double wgt=detector_weights->at(*it);
				G3TimestreamConstPtr ts = timestreams->at(*it);
				bin_tod_pol(*ts,0,0,pnt,wgt,*T,*T,*T);
				if (is_weight_){
					bin_w(*ts,pnt,wgt,*W);
				}
			}
		} else {
			if (!include_pol_rotation_){
				for (auto it=ts_ids.begin(); it!=ts_ids.end();it++){
					if (timestreams->find(*it) == timestreams->end()){
						log_info("Skipping %s", it->c_str());
						continue;
					}

					if (det_pnting->find(*it) == det_pnting->end()){
						log_fatal("No pointing information for %s",
							  it->c_str());
					}

					const std::vector<int>&pnt=det_pnting->at(*it);
					G3MapDouble::const_iterator wgtit = detector_weights->find(*it);
					double wgt = (wgtit == detector_weights->end()) ? 0 : wgtit->second;
					G3TimestreamConstPtr ts = timestreams->at(*it);
					const BolometerProperties & bp=bolo_props_->at(*it);
					bin_tod_pol(*ts,bp.pol_angle,bp.pol_efficiency,
						    pnt,wgt,*T,*Q,*U);
					if (is_weight_){
						bin_w_pol(*ts,bp.pol_angle,bp.pol_efficiency,
							  pnt,wgt,*W);
					}
				}
			}else{
				G3VectorQuatConstPtr trans=frame->Get<G3VectorQuat>(trans_key_);
				for (auto it=ts_ids.begin(); it!=ts_ids.end();it++){
					if (timestreams->find(*it) == timestreams->end()){
						log_info("Skipping %s", it->c_str());
						continue;
					}
					if (det_pnting->find(*it) == det_pnting->end()){
						log_fatal("No pointing information for %s",
							  it->c_str());
					}

					const std::vector<int>&pnt=det_pnting->at(*it);
					double wgt=detector_weights->at(*it);
					G3TimestreamConstPtr ts = timestreams->at(*it);
					const BolometerProperties & bp=bolo_props_->at(*it);
					std::vector<double> rots;
					get_detector_rotation( bp.x_offset, bp.y_offset,
							       *trans, rots);
					bin_tod_pol_rot(*ts,bp.pol_angle,bp.pol_efficiency,
							pnt,wgt,rots,*T,*Q,*U);
					if (is_weight_){
						bin_w_pol_rot(*ts,bp.pol_angle,bp.pol_efficiency,
							      pnt,wgt,rots,*W);
					}
				}
			}
		}
		out.push_back(frame);     
	} else if(frame->type == G3Frame::Map){
		g3_assert(frame->Has("Id"));
		if (map_id_ == frame->Get<G3String>("Id")->value){
			std::vector<std::string> maps_to_bin = {"T", "Q", "U", "Wpol", "Wunpol"};
			for (auto s = maps_to_bin.begin(); s < maps_to_bin.end(); s++){
				if (! frame->Has(*s) ){
					continue;
				}
				if (*s == "Wpol" || *s == "Wunpol"){
					W = boost::make_shared<G3SkyMapWeights>(
						*(frame->Get<G3SkyMapWeights>(*s)));
				} else {
					if( *s == "T") {
						T = frame->Get< G3SkyMap >(*s)->Clone(0);
						T->units = G3Timestream::None;
					} else if (*s == "Q") {
						Q = frame->Get< G3SkyMap >(*s)->Clone(0);
						Q->units = G3Timestream::None;
					} else if (*s == "U") {
						U = frame->Get< G3SkyMap >(*s)->Clone(0);
						U->units = G3Timestream::None;
					}
				}
				frame->Delete(*s);//so when we add it back later it doesn't exist
			}
			g3_assert( T );
			if (W) {
				is_weight_ = true;
			} else {
				is_weight_ = false;
			}
			if (Q){
				g3_assert(U);
				is_pol_ = true;
			} else {
				is_pol_ = false;
			}
			map_frame_ = frame;
		} else{
			out.push_back(frame);
		}
	} else if (frame->type == G3Frame::Calibration && frame->Has(bolo_props_key_)){
		bolo_props_ = frame->Get<BolometerPropertiesMap>(bolo_props_key_);
		out.push_back(frame);
	} else if (frame->type != G3Frame::EndProcessing){
		out.push_back(frame);
	}

	//Handle map emmission
	if ( map_frame_ && ( frame->type == G3Frame::EndProcessing  || 
			     (frame->type == G3Frame::Scan && frame->Has("EmitMap") && 
			      frame->Get<G3Int>("EmitMap")))){
		if(T){G3FrameObjectPtr x = boost::dynamic_pointer_cast<G3FrameObject>(T); map_frame_->Put("T", x);}
		if(Q){G3FrameObjectPtr x = boost::dynamic_pointer_cast<G3FrameObject>(Q); map_frame_->Put("Q", x);}
		if(U){G3FrameObjectPtr x = boost::dynamic_pointer_cast<G3FrameObject>(U); map_frame_->Put("U", x);}
		if(W){
			G3FrameObjectPtr x =
			    boost::dynamic_pointer_cast<G3FrameObject>(W);
			if (is_pol_){map_frame_->Put("Wpol", x);}
			else {map_frame_->Put("Wunpol", x);}
		}
		out.push_back(map_frame_);
		map_frame_ = NULL;
	}
	if (frame->type == G3Frame::EndProcessing){
		out.push_back(frame);
	}
}


void SimulatedTimestreamFiller::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if (frame->type == G3Frame::Map && frame->Get<G3String>("Id")->value == map_id_){
		if (frame->Has("T")){
			T = frame->Get<G3SkyMap>("T");
		}
		if (frame->Has("Q")){
			Q = frame->Get<G3SkyMap>("Q");
		}
		if (frame->Has("U")){
			U = frame->Get<G3SkyMap>("U");
		}
		g3_assert(T);
		is_pol_ = false;
		if (Q){
			is_pol_=true;
			g3_assert(U);
		}
	} else if (frame->type == G3Frame::Calibration){
		if (frame->Has(bolo_props_key_)){
			bolo_props_ = frame->Get<BolometerPropertiesMap>(bolo_props_key_);
		}
		out.push_back(frame);
	} else if (frame->type == G3Frame::Scan){
		g3_assert(bolo_props_);
		g3_assert(T);
		G3MapVectorIntConstPtr pnt = frame->Get<G3MapVectorInt>(pointing_key_);
		G3TimestreamConstPtr ts_with_times = 
			frame->Get<G3Timestream>(ts_to_get_sample_rate_key_);

		//create a list of ids to use for map making
		std::vector<std::string> ts_ids;
		if (ts_lst_key_ != "" && frame->Has(ts_lst_key_) ) {
			auto ts_id_lst = frame->Get<G3VectorString>(ts_lst_key_);
			ts_ids = std::vector<std::string>(ts_id_lst->begin(),ts_id_lst->end());
		} else {
			for (auto it=pnt->begin(); it!=pnt->end(); it++){
				ts_ids.push_back(it->first);
			}
		}


		//create our timestreams
		G3TimestreamMapPtr out_ts(new G3TimestreamMap());
		for (auto it=ts_ids.begin(); it != ts_ids.end(); it++){
			(*out_ts)[*it] = G3TimestreamPtr(
				new G3Timestream(ts_with_times->size(), 0));
		}

		if (! include_pol_rot_){
			if (is_pol_) {
				for (auto it=ts_ids.begin(); it != ts_ids.end(); it++){
					const BolometerProperties & bp = 
						bolo_props_->at((*it));
					fill_tod_pol(pnt->at(*it),*T,*Q,*U,
						     bp.pol_angle, bp.pol_efficiency,
						     *((*out_ts)[(*it)]) );
				}
			} else {
				for (auto it=ts_ids.begin(); it != ts_ids.end(); it++){
					const BolometerProperties & bp = 
						bolo_props_->at((*it));
					fill_tod_pol(pnt->at(*it),*T,*T,*T,
						     0, 0,
						     *((*out_ts)[(*it)]) );
				}
			}
		} else {
			if (is_pol_) {
				G3VectorQuatConstPtr trans=frame->Get<G3VectorQuat>(trans_key_);
				for (auto it=ts_ids.begin(); it != ts_ids.end(); it++){
					const BolometerProperties & bp = 
						bolo_props_->at((*it));
					std::vector<double> rots;
					get_detector_rotation(bp.x_offset, bp.y_offset,
							      *trans, rots);
					fill_tod_pol_rot(pnt->at(*it),*T,*Q,*U,
							 rots,
							 bp.pol_angle, bp.pol_efficiency,
							 *((*out_ts)[(*it)]) );
				}
			} else {
				log_fatal("You can't rotate nothing");
			}
		}
		for (auto it=out_ts->begin(); it != out_ts->end(); it++){
			it->second->start = ts_with_times->start;
			it->second->stop = ts_with_times->stop;
		}

		frame->Put( out_ts_key_, out_ts);

		out.push_back(frame);
	} else {
		out.push_back(frame);
	}
}




void SkyMapInterpTodFiller::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if (frame->type == G3Frame::Map && frame->Get<G3String>("Id")->value == map_id_){
		if (frame->Has("T")){
			T = frame->Get<G3SkyMap>("T");
		}
		if (frame->Has("Q")){
			Q = frame->Get<G3SkyMap>("Q");
		}
		if (frame->Has("U")){
			U = frame->Get<G3SkyMap>("U");
		}
		g3_assert(T);
		if (Q){
			g3_assert(U);
		} else {
			log_fatal("If you are running interpolated sims we only "
				  "support polarized sims");
		}
	} else if (frame->type == G3Frame::Calibration){
		if (frame->Has(bolo_props_key_)){
			bolo_props_ = frame->Get<BolometerPropertiesMap>(bolo_props_key_);
		}
		out.push_back(frame);
	} else if (frame->type == G3Frame::Scan){
		g3_assert(bolo_props_);
		g3_assert(T);
		G3VectorStringConstPtr ids;
		if (frame->Has(valid_ids_)){
			ids = frame->Get<G3VectorString>(valid_ids_);
		}else {
			G3VectorStringPtr id_tmp(new G3VectorString());
			for (auto it=bolo_props_->begin(); it!=bolo_props_->end(); it++){
				id_tmp->push_back(it->first);
			}
			ids = id_tmp;
		} 

		G3TimestreamConstPtr ts_with_times = 
			frame->Get<G3Timestream>(ts_to_get_sample_rate_key_);
		size_t ts_len = ts_with_times->size();
		G3TimestreamMapPtr out_ts(new G3TimestreamMap());
		for (auto it=ids->begin(); it!=ids->end(); it++){
			(*out_ts)[*it] = 
				G3TimestreamPtr(new G3Timestream(ts_len, 0));
		}

		G3MapVectorDoubleConstPtr detector_alpha;
		G3MapVectorDoubleConstPtr detector_delta;
		G3VectorQuatConstPtr trans;

		// checks if we have already calculated the detector pointing
		bool need_to_calc_pnt = true;
		if ( pre_computed_pointing_ &&
		     frame->Has(detector_alpha_key_) && frame->Has(detector_delta_key_) ){
			detector_alpha = frame->Get<G3MapVectorDouble>(detector_alpha_key_);
			detector_delta = frame->Get<G3MapVectorDouble>(detector_delta_key_);
			need_to_calc_pnt = false;
		} else {
			trans = frame->Get<G3VectorQuat>(trans_key_);
		}

		for (size_t ts_ind=0; ts_ind < ids->size(); ts_ind++){
			std::string id = ids->at(ts_ind);
			const BolometerProperties & bp = bolo_props_->at((id));
			std::vector<double> alpha(trans->size());
			std::vector<double> delta(trans->size());
			std::vector<double> rots(ts_len, 0);
			if (include_pol_rot_) {
				get_detector_rotation(bp.x_offset, bp.y_offset,
						      *trans, rots);
			}
			if (! need_to_calc_pnt){
				alpha = std::vector<double>(detector_alpha->at((id)));
				delta = std::vector<double>(detector_delta->at((id)));
			} else {
				get_detector_pointing(bp.x_offset, bp.y_offset,
						      *trans, T->coord_ref, alpha, delta);
			}

			fill_tod_pol_rot_interp_2d(alpha, delta, *T,*Q,*U,
						   rots,
						   bp.pol_angle, bp.pol_efficiency,
						   *((*out_ts)[(id)]) );
			
		}
		for (auto it=out_ts->begin(); it != out_ts->end(); it++){
			it->second->start = ts_with_times->start;
			it->second->stop = ts_with_times->stop;
		}

		frame->Put( out_ts_key_, out_ts);
		out.push_back(frame);
	} else {
		out.push_back(frame);
	}
}


namespace bp = boost::python;

EXPORT_G3MODULE("mapmaker", MapPointingCalculator,
		(bp::init<std::string, std::string, std::string, std::string,
		 std::string, std::string, std::string> (
			 (
				 bp::arg("map_id"), 
				 bp::arg("ts_map_key"),
				 bp::arg("bolo_props_key"),
				 bp::arg("trans_key"),
				 bp::arg("detector_alpha_key"),
				 bp::arg("detector_delta_key"),
				 bp::arg("pointing_store_key")))), 
		MAP_POINTING_CALCULATOR_DOCSTR);

EXPORT_G3MODULE("mapmaker", MapBinner,
		(bp::init<std::string, std::string, std::string, std::string,
		 std::vector<std::string>, std::string, std::string, bool> (
			 (
				 bp::arg("map_id"), 
				 bp::arg("ts_map_key"),
				 bp::arg("pointing_key"),
				 bp::arg("weight_key") = "",
				 bp::arg("bolo_ids_to_map") = std::vector<std::string>(),
				 bp::arg("bolo_props_key") = "BolometerProperties",
				 bp::arg("trans_key") = "RaDecTransform",
				 bp::arg("include_pol_rotation")=false
				 ))
			), 
		MAP_BINNER_DOCSTR);

PYBINDINGS("mapmaker"){
	bp::class_<SimulatedTimestreamFiller,
		bp::bases<G3Module>,
		boost::shared_ptr<SimulatedTimestreamFiller>,
		boost::noncopyable
		>("SimulatedTimestreamFiller", bp::init<std::string, std::string, std::string, 
		  std::string, std::string, std::string, bool, std::string>
		  ("Doc", 
		   (bp::arg("map_id"),
		    bp::arg("sim_pointing_key"),
		    bp::arg("bolo_props_key"),
		    bp::arg("ts_to_get_sample_rate_key"),
		    bp::arg("trans_key"),
		    bp::arg("ts_lst_key"),
		    bp::arg("include_pol_rot"),
		    bp::arg("out_ts_key"))))
		  .def("Process", &SimulatedTimestreamFiller::Process);


	bp::class_<SkyMapInterpTodFiller,
		bp::bases<G3Module>,
		boost::shared_ptr<SkyMapInterpTodFiller>,
		boost::noncopyable
		>("SkyMapInterpTodFiller", bp::init<std::string, std::string, std::string,
		  std::string, bool,std::string, std::string, std::string, std::string>
		  ("Doc", 
		   (bp::arg("map_id"),
		    bp::arg("bolo_props_key"),
		    bp::arg("ts_to_get_sample_rate_key"),
		    bp::arg("trans_key"),
		    bp::arg("include_pol_rot"),
		    bp::arg("valid_ids"),
		    bp::arg("detector_alpha_key"),
		    bp::arg("detector_delta_key"),
		    bp::arg("out_ts_key"))))
		  .def("Process", &SkyMapInterpTodFiller::Process);

}
