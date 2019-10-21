#include <pybindings.h>
#include <timestreamflagging/removeflags.h>


void RemoveFlagged::Process(G3FramePtr frame, std::deque<G3FramePtr> &out){
	if (frame->type == G3Frame::Scan){
                g3_assert(frame->Has(in_key_));
                g3_assert(frame->Has(flag_key_));
                G3TimestreamMapConstPtr tsin = frame->Get<G3TimestreamMap>(in_key_);
                G3TimestreamMapPtr tsout = boost::make_shared<G3TimestreamMap>( *tsin );
                G3MapVectorStringConstPtr flagged_channels = 
			frame->Get<G3MapVectorString>(flag_key_);
                for (auto k =flagged_channels->begin(); k != flagged_channels->end(); k++){
			auto pos = tsout->find(k->first);
                        if (pos != tsout->end()){
                                tsout->erase(pos);
			}
                }
                frame->Put( out_key_, tsout);
        }
        out.push_back(frame);
}

EXPORT_G3MODULE("timestreamflagging", RemoveFlagged,
		(init<std::string, std::string,  std::string>(
		   (arg("input_ts_key"), arg("input_flag_key"),arg("output_ts_key")))), 
		"");
