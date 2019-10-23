import rawdumptools
import sys

path = sys.argv[1]

# preemptively look for rawdumps to process in the data directory
while True:
    data_files = os.listdir(path)
    filename = [file for file in data_files if '.lock' not in file][0]
    lock_filename = filename + '.lock'

    # change name to '.lock' to prevent other worker processes from reading
    shutil.move(filename, lock_filename)

    # run the processing
    rawdumptools.process(path + lock_filename, path + '/processed/' + filename + '.pkl')

