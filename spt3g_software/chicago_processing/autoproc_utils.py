import os, sys
import shutil
import subprocess
import tempfile
import pandas
from spt3g.cluster.condor_tools import condor_submit_baked

# Utilities for creating and launching  autoprocessingjobs, as well as
# introspecting the autoprocessing dependencies

softwareroot = os.getenv('SPT3G_SOFTWARE_PATH')


def get_autoproc_args(as_dict=False, **kwargs):
    '''
    Check that all arguments required for running autoprocessing jobs are set.

    Arguments
    ---------
    as_dict : bool
        If True, return the arguments below as a dctionary rather than as a tuple
    calresults : str
        Path to autoprocessing output directory
    proxycert : str
        Path to the grid proxy certificate, required for grid file transfers
    atpole : bool
        If True, this instance is running at pole, otherwise on the northern grid
    toolset : str
        Version of clustertools to link against.

    Returns
    -------
    calresults, proxycert, atpole, toolset
        As extracted from the input kwargs dictionary
    kwargs : dict
        Remaining arguments to be passed on elsewhere
    '''
    calresults = kwargs.pop('calresults', None)
    proxycert = kwargs.pop('proxycert', None)
    atpole = kwargs.pop('atpole', None)
    toolset = kwargs.pop('toolset', None)
    for arg in [calresults, atpole, toolset]:
        if arg is None:
            raise RuntimeError(
                'One or more required arguments are missing for running '
                'autoprocessing jobs. Did you run set_autoproc_args() ?'
            )
    if not os.path.exists(calresults):
        raise OSError("Missing calresults directory {}".format(calresults))
    if not atpole and (not proxycert or not os.path.exists(proxycert)):
        raise OSError("Missing grid proxy certificate {}".format(proxycert))

    if as_dict:
        ret = {
            'calresults': calresults,
            'proxycert': proxycert,
            'atpole': atpole,
            'toolset': toolset,
        }
        return ret, kwargs
    return calresults, proxycert, atpole, toolset, kwargs


# General-purpose processing execution functions: condor_submit() and
# execute_locally(). They both have identical signatures, assume that the
# script has a signature like script.py -o output <input files>, and take
# keyword arguments:
# - memory: expected memory consumption in MB
# - outfile: path to the output file under calresults/calibration, if different
#   from source/obsid.g3


def condor_submit(source, obsid, script, infiles, *args, **kwargs):
    '''Submit a condor job'''
    calresults, proxycert, atpole, toolset, kwargs = get_autoproc_args(**kwargs)
    script = os.path.join(softwareroot, script)
    logdir = os.path.join(calresults, 'logs', source, str(obsid))
    if 'outfile' in kwargs:
        outfile = os.path.join(calresults, 'calibration', kwargs['outfile'])
    else:
        outfile = os.path.join(calresults, 'calibration', source, '%d.g3' % obsid)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    submit = [
        'executable = autoproc_job.sh',
        'universe = vanilla',
        'error = %s' % os.path.join(logdir, 'error.log'),
        'output = %s' % os.path.join(logdir, 'out.log'),
        'log = %s' % os.path.join(logdir, 'condor.log'),
        '+WANT_RCC_ciconnect = True',
        '+ProjectName = "spt.all"',
        'request_cpus = 1',
        'transfer_output_files = ""',
        'should_transfer_files = YES',
        'run_as_owner = %s' % ('True' if atpole else 'False'),
        'when_to_transfer_output = ON_EXIT',
        'transfer_input_files = %s,%s' % (script, proxycert),
        'transfer_executable = True',
        'arguments = "\'%s\' \'%s\' \'%s\' \'%s\'"'
        % (os.path.basename(script), outfile, toolset, ' '.join(infiles)),
    ]

    # Add proxy cert if we have one
    if proxycert:
        submit.append('x509userproxy = %s' % proxycert)

    # Construct environment variables
    env = 'SPT_SCRIPT_ARGS=\'%s\'' % ' '.join(args)
    if atpole:
        submit.append('priority = 100')
        submit.append('getenv = True')
        env += ' SPT_NO_CVMFS=1'
    submit.append('environment = "%s"' % env)

    # Estimate resource consumption if we can
    if 'memory' in kwargs:
        submit.append('request_memory = %d' % kwargs['memory'])

    # Get the kinds of machines we want when on the grid
    if not atpole:
        requirements = '(OSGVO_OS_STRING == "RHEL 6" || OSGVO_OS_STRING == "RHEL 7")'
        requirements += '&& (HAS_CVMFS_spt_opensciencegrid_org)'

        # Site blacklist as a result of bad CVMFS configuration
        requirements += '&& (GLIDEIN_ResourceName =!= "SPRACE")'  # Top two are Brazil
        requirements += '&& (GLIDEIN_ResourceName =!= "GridUNESP_CENTRAL")'
        requirements += '&& (GLIDEIN_ResourceName =!= "NPX")'

        if 'memory' in kwargs and kwargs['memory'] > 16384:
            requirements += '|| SPT_BIGMEM =?= True'
        submit.append('Requirements = (%s)' % requirements)

    disk = 100 * 1024 * 1024  # 100 MB of scratch at minimum
    for f in infiles:  # Plus all the input files
        try:
            disk += os.stat(f).st_size
        except Exception as e:
            print(
                "Error submitting job %s/%d: getting file size: %s"
                % (source, obsid, str(e)),
                file=sys.stderr,
            )
            return None
    submit.append('request_disk = %d' % (disk / 1024))

    submit.append('queue')
    try:
        jobid, result = condor_submit_baked(submit)
        if jobid is None:
            print(
                'Error submitting job %s/%d: no jobid: %s' % (source, obsid, result),
                file=sys.stderr,
            )
        return jobid
    except subprocess.CalledProcessError as e:
        print(
            'Error submitting job %s/%d: condor error: %s'
            % (source, obsid, e.output.decode() if e.output else None),
            file=sys.stderr,
        )
        return None


def execute_locally(source, obsid, script, infiles, *args, **kwargs):
    '''
    Local execution for processing that is so simple that the overhead of Condor
    makes no sense
    '''
    calresults, proxycert, atpole, toolset, kwargs = get_autoproc_args(**kwargs)
    script = os.path.join(softwareroot, script)
    if 'outfile' in kwargs:
        out = os.path.join(calresults, 'calibration', kwargs['outfile'])
    else:
        out = os.path.join(calresults, 'calibration', source, '%d.g3' % obsid)

    tmp = tempfile.mkstemp(suffix='.g3')
    tmpout = tmp[1]
    command = [os.path.join(softwareroot, script), '-o', tmpout]
    command += args
    command += infiles
    try:
        result = subprocess.check_output(command)
        # Use GridFTP to copy back if there is a proxy cert defined
        if proxycert:
            copyenv = os.environ.copy()
            copyenv['X509_USER_PROXY'] = proxycert
            result = subprocess.check_output(
                [
                    'globus-url-copy',
                    '-cd',
                    'file://' + tmpout,
                    'gsiftp://gridftp.grid.uchicago.edu:2811/cephfs' + out,
                ],
                env=copyenv,
            )
        else:  # The local case (usually at pole)
            if not os.path.exists(os.path.dirname(out)):
                os.makedirs(os.path.dirname(out))
            shutil.copy(tmpout, out)
            os.chmod(out, 0o644)
    except subprocess.CalledProcessError as e:
        print("Error running job %s/%d" % (source, obsid), file=sys.stderr)
        print(
            'Error during local script execution (%s): %s'
            % (' '.join(command), e.output.decode() if e.output else None),
            file=sys.stderr,
        )
        return 'failed'
    finally:
        os.unlink(tmpout)
        os.close(tmp[0])

    return 'complete'


def make_calframe(src, obs, files, **kwargs):
    '''Assemble outputs from several jobs into a single calframe'''
    args, _ = get_autoproc_args(as_dict=True, **kwargs)
    return execute_locally(
        src.split('/')[0],
        obs,
        'calibration/scripts/build_cal_frame.py',
        files,
        '-c',
        os.path.join(args['calresults'], 'config'),
        '-t',
        (
            pandas.Timestamp('2017-01-01 00:00:00', tz='UTC')
            + pandas.Timedelta(seconds=obs)
        ).isoformat(),
        outfile='calframe/%s/%d.g3' % (src.split('/')[0], obs),
        **args,
    )


def get_script(a, b, script, files, *args, return_args=False, **kwargs):
    '''
    A dummy function that takes the place of the above functions.
    This is used when performing introspection to get the name of the script
    that would be called, and its arguments.

    Most of the arguments are unused in this function, and only exist
    to maintain the API with the function that runs the scripts.
    '''
    if not return_args:
        return script
    cmd = [os.path.join(softwareroot, script), '-o', '<output-file>']
    cmd += [str(x) for x in args]
    cmd += ['<input-files>']
    return ' '.join(cmd)


def get_calframe_script(a, b, script, return_args=False, **kwargs):
    '''
    Another dummy function that just spits out information about `make_calframe`
    '''
    if return_args:
        return 'make_calframe (defined in {})'.format(os.path.abspath(__file__))
    return 'make_calframe'


def make_obs_type_dict():
    '''
    Classify sources by their "type"- i.e. how they are used in analysis.
    '''
    return {
        'fluxandpointing': ['RCW38-pixelraster', 'MAT5A-pixelraster'],
        'opacity': ['RCW38', 'MAT5A'],
        'planets': ['mars', 'saturn', 'moon'],  # "planets"
    }


# Processing information: maps source name to processing function and
# (optional) dependencies list
def make_scripts_dict(for_dependency_graph=False, **kwargs):
    '''
    Returns a dictionary which maps autoprocessing targets to the
    scripts that construct them, and their dependencies.
    If `for_dependency_graph` is False, all other arguments are required.

    Set get_autoproc_args() for more details.
    '''
    if not for_dependency_graph:
        args, _ = get_autoproc_args(as_dict=True, **kwargs)
        global condor_submit, execute_locally, make_calframe
    else:
        condor_submit = get_script
        execute_locally = get_script
        make_calframe = get_calframe_script

    make_calframe_args = lambda src, obs, files: make_calframe(
        src, obs, files, **args
    )

    scripts = {
        'calibrator': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/analyze_calibrator.py',
                files,
                memory=8192,
                **args,
            ),
            ['elnod', 'boloproperties', 'tod'],
        ),
        'elnod': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/analyze_elnod.py',
                files,
                memory=4096,
                **args,
            ),
            ['tod'],
        ),
        # RCW38 fast points
        'RCW38-pixelraster/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makemaps.py',
                files,
                memory=30720,
                **args,
            ),
            ['calibrator', 'elnod', 'tod', 'boloproperties'],
        ),
        'RCW38-pixelraster': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/fit_fluxandpointing.py',
                files,
                memory=30720,
                **args,
            ),
            ['calibrator', 'boloproperties', 'RCW38-pixelraster/maps'],
        ),
        # RCW38 very fast points
        'RCW38/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-k',
                memory=4096,
                **args,
            ),
            ['elnod', 'calibrator', 'boloproperties', 'RCW38-pixelraster', 'tod'],
        ),
        'RCW38': (
            lambda src, obs, files: execute_locally(
                src, obs, 'calibration/scripts/fluxpointcal/opacity.py', files, **args
            ),
            ['boloproperties', 'RCW38-pixelraster', 'RCW38/maps'],
        ),
        # MAT5A fast points
        'MAT5A-pixelraster/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makemaps.py',
                files,
                memory=30720,
                **args,
            ),
            ['calibrator', 'elnod', 'tod', 'boloproperties'],
        ),
        'MAT5A-pixelraster': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/fit_fluxandpointing.py',
                files,
                '-s',
                'MAT5A',
                memory=30720,
                **args,
            ),
            ['calibrator', 'boloproperties', 'MAT5A-pixelraster/maps'],
        ),
        # MAT5A very fast points
        'MAT5A/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-k',
                memory=4096,
                **args,
            ),
            ['elnod', 'calibrator', 'boloproperties', 'MAT5A-pixelraster', 'tod'],
        ),
        'MAT5A': (
            lambda src, obs, files: execute_locally(
                src,
                obs,
                'calibration/scripts/fluxpointcal/opacity.py',
                files,
                '-s',
                'MAT5A',
                **args,
            ),
            ['boloproperties', 'MAT5A-pixelraster', 'MAT5A/maps'],
        ),
        # Other source observations
        'mars-pixelraster/singlebolomaps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makemaps.py',
                files,
                memory=30720,
                **args,
            ),
            ['elnod', 'calibrator', 'tod', 'boloproperties'],
        ),
        'saturn-pixelraster/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makemaps.py',
                files,
                memory=30720,
                **args,
            ),
            ['elnod', 'calibrator', 'tod', 'boloproperties'],
        ),
        'saturn/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-r',
                '0.2',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        'saturn': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/fit_fluxandpointing.py',
                files,
                memory=30720,
                **args,
            ),
            ['calibrator', 'saturn/maps'],
        ),
        # Observations that do not feed into calframe generation and can depend
        # on final calframe outputs
        'CenA-pixelraster/singlebolomaps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makemaps.py',
                files,
                '-t',
                '-p',
                '-r',
                '0.25',
                '-x',
                '0.5',
                '-y',
                '0.5',
                memory=8192,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        'CenA-pixelraster/coaddmaps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-p',
                '-r',
                '0.25',
                '-x',
                '0.5',
                '-y',
                '0.5',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        'mars-pixelraster/maps': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-r',
                '0.2',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        '0537-441': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-r',
                '0.2',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        '0537-441-pixelraster': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-r',
                '0.2',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        'PMNJ0210-5101': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-r',
                '0.2',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        'PMNJ0210-5101-pixelraster': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/fluxpointcal/makecoadd.py',
                files,
                '-r',
                '0.2',
                memory=4096,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        'noise': (
            lambda src, obs, files: condor_submit(
                src,
                obs,
                'calibration/scripts/analyze_noise.py',
                files,
                '--cal-sn-threshold',
                '10',
                memory=40000,
                **args,
            ),
            ['calframe', 'tod'],
        ),
        # Meta-products
        'calframe': (
            make_calframe_args,
            ['elnod', 'calibrator', 'boloproperties', 'opacity'],
        ),
        # calframes for parts of the calibration chain
        'elnod/calframe': (make_calframe_args, ['elnod', 'boloproperties']),
        'calibrator/calframe': (
            make_calframe_args,
            ['elnod', 'calibrator', 'boloproperties'],
        ),
        'RCW38-pixelraster/calframe': (
            make_calframe_args,
            ['elnod', 'calibrator', 'boloproperties', 'RCW38-pixelraster'],
        ),
        'MAT5A-pixelraster/calframe': (
            make_calframe_args,
            ['elnod', 'calibrator', 'boloproperties', 'MAT5A-pixelraster'],
        ),
        # NB: 'boloproperties' is inserted into the DB by hand
        'default': None,  # (lambda src, obs, files: None, ['tod', 'calibrator', 'RCW38'])
    }
    return scripts


##########################################################################
# Code below is for parsing the dependency graph of autoprocessing


def reduce_similar_obs(obsname):
    '''
    Certain sources are treated identically in autoprocessing, so we replace
    them with a single name to reduce the complexity of the dependency graph.
    '''
    similarobs = {
        'RCW38': '[flux_cal_src]',
        'MAT5A': '[flux_cal_src]',
        '0537-441': '[focus_quasar]',
        'PMNJ0210-5101': '[focus_quasar]',
    }
    source = obsname.replace('-pixelraster', '')
    source = source.split('/')[0]
    if source in similarobs:
        obsname = obsname.replace(source, similarobs[source])
    return obsname


def add_script_name(scripts, graph):
    '''
    Add script names to the graph nodes.
    '''
    for node in graph.nodes():
        node = graph.get_node(node)
        node.attr['label'] = '<<b>{}</b>>'.format(node.name)
    for target in scripts:
        targname = reduce_similar_obs(target)
        if not graph.has_node(targname):
            continue
        node = graph.get_node(targname)
        script_name = os.path.basename(scripts[target][0](None, None, None))
        node.attr['label'] = '<<b> {} </b><br/><i>{}</i>>'.format(
            node.name, script_name
        )
    return graph


def build_autoproc_depdict(scripts):
    '''
    Build a dictionary parseable by pygraphviz that describes the
    autoprocessing dependencies (excluding planets).
    '''
    depdict = {}
    planets = make_obs_type_dict()['planets']
    for target in scripts:
        if scripts[target] is None:
            continue
        if any([planet in target for planet in planets]):
            # drop if this is a planet job
            continue
        if 'calframe' in target:
            # We don't care about the different kinds of calframe
            target = 'calframe'
        targname = reduce_similar_obs(target)
        if 'planet' in targname:
            continue
        for dep in scripts[target][1]:
            if 'calframe' in dep:
                dep = 'calframe'
            depname = reduce_similar_obs(dep)
            if depname == 'opacity':
                # Opacity as a dependency is the same a flux_cal_src as a target
                depname = '[flux_cal_src]'
            if depname not in depdict:
                depdict[depname] = {targname: None}
            else:
                depdict[depname][targname] = None
    depdict.pop('tod')  # Everything depends on TOD, and that is obvious
    return depdict


def build_planet_depdict(scripts):
    '''
    Build a dictionary describing the dependencies of the planet
    processing , without the full autoprocessing dependencies.
    '''
    depdict = {}
    planets = make_obs_type_dict()['planets']
    for target in scripts:
        if scripts[target] is None:
            continue
        if not any([planet in target for planet in planets]):
            # drop if this isn't a planet job
            continue
        targname = target
        for dep in scripts[target][1]:
            depname = reduce_similar_obs(dep)
            if depname == 'opacity':
                depname = '[flux_cal_src]'
            if depname not in depdict:
                depdict[depname] = {targname: None}
            else:
                depdict[depname][targname] = None
    depdict.pop('tod')
    return depdict


def build_graph(scripts, depdict, dotfile=None, imgfile=None):
    '''
    Convert a dependency dictionary to a graph.
    If `dotfile` is not None, write the graph as a text file to `dotfile`.
    If `imgfile` is not None, write the graph as an image to `imgfile`.
    '''
    try:
        import pygraphviz as pgv
    except ImportError as err:
        print('pygraphviz can be acquired by running `pip install --user pygraphviz`')
        raise err
    graph = pgv.AGraph(depdict, strict=True, directed=True)
    graph = add_script_name(scripts, graph)
    if imgfile is not None:
        graph.layout(prog='dot')
        graph.draw(imgfile)
    if dotfile is not None:
        graph.write(dotfile)
    return graph


def make_default_dependency_graphs():
    '''
    Spit out some nice pretty pictures for use in documentation
    '''
    # TODO: assign reasonable locations for these things!
    # But, we can punt on that until we actually integrate this with the docs.
    scripts = make_scripts_dict(for_dependency_graph=True)
    gr = build_graph(
        scripts,
        build_autoproc_depdict(scripts),
        imgfile='autograph.png',
        dotfile='autograph.txt',
    )
    build_graph(scripts, build_planet_depdict(scripts), imgfile='planetgraph.png')


if __name__ == '__main__':
    make_default_dependency_graphs()
