import os

def get_test_files_path():
    assert('SPT3G_SOFTWARE_BUILD_PATH' in os.environ)
    return os.environ['SPT3G_SOFTWARE_BUILD_PATH']+'/testdata/mapmaker_test_data/'
