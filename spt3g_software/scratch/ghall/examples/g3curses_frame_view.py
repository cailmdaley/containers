from spt3g import core
from spt3g.util import G3Viewer

v1 = core.G3VectorString(["foo", "bar", "baz"])
v2 = core.G3VectorString(["some", "more", "values"])
v3 = core.G3VectorDouble([3.14, 2.72, 144])
v4 = core.G3VectorDouble([1.0, 4.0, 9.0])

m1 = core.G3MapVectorString()
m1["first string vect"] = v1
m1["second string vect"] = v2

m2 = core.G3MapVectorDouble()
m2["first double vect"] = v3
m2["second double vect"] = v4

f1 = core.G3Frame(core.G3FrameType.Map)
ts1 = core.G3Timestream([1,2,3,4])
f1['A G3Timestream'] = ts1
f1['A map of string vectors'] = m1

f2 = core.G3Frame(core.G3FrameType.Scan)
ts2 = core.G3Timestream([5,6,7,8])
f2['A Different G3Timestream'] = ts2
f2['A map of double vectors'] = m2

d = [f1, f2]

v = G3Viewer(d, current_year=1993)
