#pragma once
#include <G3Module.h>
#include <G3Logging.h>
#include <G3Timestream.h>
#include <G3Map.h>

class RemoveFlagged : public G3Module{
public:
        RemoveFlagged(std::string in_ts_map_key, 
		      std::string flag_key,
		      std::string out_ts_map_key
		) :
		in_key_(in_ts_map_key), out_key_(out_ts_map_key),
		flag_key_(flag_key){}
        void Process(G3FramePtr frame, std::deque<G3FramePtr> &out);
private:
	std::string in_key_;
	std::string out_key_;
	std::string flag_key_;
        SET_LOGGER("RemoveFlagged");
};
G3_POINTERS(RemoveFlagged);


