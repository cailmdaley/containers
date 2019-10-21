#!/bin/bash

#curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=NLHARR@berkeley.edu&password=asdfasdfasdfasdf'
#curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2016-06-27--2016-06-29/format/3le' > 3le1606
#curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org//ajaxauth/logout


curl -c cookies.txt -b cookies.txt https://www.space-track.org/ajaxauth/login -d 'identity=NLHARR@berkeley.edu&password=asdfasdfasdfasdf'
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-03-28--2015-03-30/format/3le' > 3le_0.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-04-20--2015-04-22/format/3le' > 3le_1.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-06-18--2015-06-20/format/3le' > 3le_2.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-07-19--2015-07-21/format/3le' > 3le_3.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-08-15--2015-08-17/format/3le' > 3le_4.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-09-02--2015-09-04/format/3le' > 3le_5.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2015-10-07--2015-10-09/format/3le' > 3le_6.tles
curl --limit-rate 100K --cookie cookies.txt 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/2016-06-27--2016-06-29/format/3le' > 3le_7.tles
curl --limit-rate 100K --cookie cookies.txt https://www.space-track.org//ajaxauth/logout

