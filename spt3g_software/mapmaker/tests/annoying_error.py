from spt3g import core, mapmaker

pipe = core.G3Pipeline()
pipe.Add(mapmaker.StokesStuffer, bolo_properties = core.G3VectorDouble())
pipe.Run()
