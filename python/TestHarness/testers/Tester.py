import platform, re, os
from TestHarness import util
from FactorySystem.MooseObject import MooseObject
from tempfile import TemporaryFile
import subprocess
from signal import SIGTERM

class Tester(MooseObject):
    """
    Base class from which all tester objects are instanced.
    """
    @staticmethod
    def validParams():
        params = MooseObject.validParams()

        # Common Options
        params.addRequiredParam('type', "The type of test of Tester to create for this test.")
        params.addParam('max_time',   300, "The maximum in seconds that the test will be allowed to run.")
        params.addParam('min_reported_time', 10, "The minimum time elapsed before a test is reported as taking to long to run.")
        params.addParam('skip',     "Provide a reason this test will be skipped.")
        params.addParam('deleted',         "Tests that only show up when using the '-e' option (Permanently skipped or not implemented).")

        params.addParam('heavy',    False, "Set to True if this test should only be run when the '--heavy' option is used.")
        params.addParam('group',       [], "A list of groups for which this test belongs.")
        params.addParam('prereq',      [], "A list of prereq tests that need to run successfully before launching this test.")
        params.addParam('skip_checks', False, "Tells the TestHarness to skip additional checks (This parameter is set automatically by the TestHarness during recovery tests)")
        params.addParam('scale_refine',    0, "The number of refinements to do when scaling")
        params.addParam('success_message', 'OK', "The successful message")

        params.addParam('cli_args',       [], "Additional arguments to be passed to the test.")
        params.addParam('allow_test_objects', False, "Allow the use of test objects by adding --allow-test-objects to the command line.")

        params.addParam('valgrind', 'NONE', "Set to (NONE, NORMAL, HEAVY) to determine which configurations where valgrind will run.")
        params.addParam('tags',      [], "A list of strings")

        # Test Filters
        params.addParam('platform',      ['ALL'], "A list of platforms for which this test will run on. ('ALL', 'DARWIN', 'LINUX', 'SL', 'LION', 'ML')")
        params.addParam('compiler',      ['ALL'], "A list of compilers for which this test is valid on. ('ALL', 'GCC', 'INTEL', 'CLANG')")
        params.addParam('petsc_version', ['ALL'], "A list of petsc versions for which this test will run on, supports normal comparison operators ('<', '>', etc...)")
        params.addParam('petsc_version_release', ['ALL'], "A test that runs against PETSc master if FALSE ('ALL', 'TRUE', 'FALSE')")
        params.addParam('slepc_version', [], "A list of slepc versions for which this test will run on, supports normal comparison operators ('<', '>', etc...)")
        params.addParam('mesh_mode',     ['ALL'], "A list of mesh modes for which this test will run ('DISTRIBUTED', 'REPLICATED')")
        params.addParam('method',        ['ALL'], "A test that runs under certain executable configurations ('ALL', 'OPT', 'DBG', 'DEVEL', 'OPROF', 'PRO')")
        params.addParam('library_mode',  ['ALL'], "A test that only runs when libraries are built under certain configurations ('ALL', 'STATIC', 'DYNAMIC')")
        params.addParam('dtk',           ['ALL'], "A test that runs only if DTK is detected ('ALL', 'TRUE', 'FALSE')")
        params.addParam('unique_ids',    ['ALL'], "A test that runs only if UNIQUE_IDs are enabled ('ALL', 'TRUE', 'FALSE')")
        params.addParam('recover',       True,    "A test that runs with '--recover' mode enabled")
        params.addParam('vtk',           ['ALL'], "A test that runs only if VTK is detected ('ALL', 'TRUE', 'FALSE')")
        params.addParam('tecplot',       ['ALL'], "A test that runs only if Tecplot is detected ('ALL', 'TRUE', 'FALSE')")
        params.addParam('dof_id_bytes',  ['ALL'], "A test that runs only if libmesh is configured --with-dof-id-bytes = a specific number, e.g. '4', '8'")
        params.addParam('petsc_debug',   ['ALL'], "{False,True} -> test only runs when PETSc is configured with --with-debugging={0,1}, otherwise test always runs.")
        params.addParam('curl',          ['ALL'], "A test that runs only if CURL is detected ('ALL', 'TRUE', 'FALSE')")
        params.addParam('tbb',           ['ALL'], "A test that runs only if TBB is available ('ALL', 'TRUE', 'FALSE')")
        params.addParam('superlu',       ['ALL'], "A test that runs only if SuperLU is available via PETSc ('ALL', 'TRUE', 'FALSE')")
        params.addParam('slepc',         ['ALL'], "A test that runs only if SLEPc is available ('ALL', 'TRUE', 'FALSE')")
        params.addParam('unique_id',     ['ALL'], "A test that runs only if libmesh is configured with --enable-unique-id ('ALL', 'TRUE', 'FALSE')")
        params.addParam('cxx11',         ['ALL'], "A test that runs only if CXX11 is available ('ALL', 'TRUE', 'FALSE')")
        params.addParam('asio',          ['ALL'], "A test that runs only if ASIO is available ('ALL', 'TRUE', 'FALSE')")
        params.addParam('depend_files',  [], "A test that only runs if all depend files exist (files listed are expected to be relative to the base directory, not the test directory")
        params.addParam('env_vars',      [], "A test that only runs if all the environment variables listed exist")
        params.addParam('should_execute', True, 'Whether or not the executable needs to be run.  Use this to chain together multiple tests based off of one executeable invocation')
        params.addParam('required_submodule', [], "A list of initialized submodules for which this test requires.")
        params.addParam('required_objects', [], "A list of required objects that are in the executable.")
        params.addParam('check_input',    False, "Check for correct input file syntax")
        params.addParam('display_required', False, "The test requires and active display for rendering (i.e., ImageDiff tests).")
        params.addParam('boost',         ['ALL'], "A test that runs only if BOOT is detected ('ALL', 'TRUE', 'FALSE')")

        # Queueing specific
        params.addParam('copy_files',         [], "Additional list of files/directories to copy when performing queueing operations")
        params.addParam('link_files',         [], "Additional list of files/directories to symlink when performing queueing operations")
        params.addParam('queue_scheduler',  True, "A test that runs only if using queue options")

        return params

    # This is what will be checked for when we look for valid testers
    IS_TESTER = True

    def __init__(self, name, params):
        MooseObject.__init__(self, name, params)
        self.specs = params
        self.outfile = None
        self.std_out = ''
        self.exit_code = 0
        self.process = None
        self.tags = params['tags']
        self.__caveats = set([])

        # Bool if test can run
        self._runnable = None

        # Initialize the status bucket class
        self.status = util.TestStatus()

        # Enumerate the buckets here so ther are easier to work with in the tester class
        self.bucket_initialized  = self.status.bucket_initialized
        self.bucket_success      = self.status.bucket_success
        self.bucket_fail         = self.status.bucket_fail
        self.bucket_diff         = self.status.bucket_diff
        self.bucket_pending      = self.status.bucket_pending
        self.bucket_finished     = self.status.bucket_finished
        self.bucket_deleted      = self.status.bucket_deleted
        self.bucket_skip         = self.status.bucket_skip
        self.bucket_silent       = self.status.bucket_silent
        self.bucket_queued       = self.status.bucket_queued
        self.bucket_waiting_processing = self.status.bucket_waiting_processing

        # Set the status message
        if self.specs['check_input']:
            self.success_message = 'SYNTAX PASS'
        else:
            self.success_message = self.specs['success_message']

        # Set up common paramaters
        self.should_execute = self.specs['should_execute']
        self.check_input = self.specs['check_input']

        if self.specs["allow_test_objects"]:
            self.specs["cli_args"].append("--allow-test-objects")

    def getTestName(self):
        """ return test name """
        return self.specs['test_name']

    def getPrereqs(self):
        """ return list of prerequisite tests this test depends on """
        return self.specs['prereq']

    def getMooseDir(self):
        """ return moose directory """
        return self.specs['moose_dir']

    def getTestDir(self):
        """ return directory this tester is located """
        return self.specs['test_dir']

    def getMinReportTime(self):
        """ return minimum time elapse before reporting a 'long running' status """
        return self.specs['min_reported_time']

    def getMaxTime(self):
        """ return maximum time elapse before reporting a 'timeout' status """
        return self.specs['max_time']

    def getRunnable(self, options):
        """ return bool and cache results, if this test can run """
        if self._runnable is None:
            self._runnable = self.checkRunnableBase(options)
        return self._runnable

    def getColor(self):
        """ return print color assigned to this tester """
        return self.status.getColor()

    def getInputFile(self):
        """ return the input file if applicable to this Tester """
        return None

    def getOutputFiles(self):
        """ return the output files if applicable to this Tester """
        return []

    def getSuccessMessage(self):
        """ return the success message assigned to this tester """
        return self.success_message

    def getStatusMessage(self):
        """ return the status message assigned to this tester """
        return self.status.getStatusMessage()

    def getStatus(self):
        """ return current enumerated tester status bucket """
        return self.status.getStatus()

    def setStatus(self, reason, bucket):
        """
        Method to set a testers status.

        Syntax:
          .setStatus('str message', <enumerated tester status bucket>)
        """
        self.status.setStatus(reason, bucket)
        return self.getStatus()

    # Method to check if a test has failed. This method will return true if a
    # tester has failed at any point during the processing of the test.
    # Note: It's possible for a tester to report false for both didFail and
    #       didPass. This will happen if the tester is in-progress for instance.
    # See didPass()
    def didFail(self):
        """
        return bool for tester failure
        see util.TestStatus for more information
        """
        return self.status.didFail()

    # Method to check for successfull test
    # Note: This method can return False until the tester has completely finished.
    #       For this reason it should be used only after the tester has completed.
    #       Instead you may want to use the didFail method which returns false
    #       only if the tester has failed at any point during the processing
    #       of that tester (e.g. after the main command has been run but before
    #       output has been tested).
    # See didFail()
    def didPass(self):
        """
        return boolean for tester successfulness
        see util.TestStatus for more information
        """
        return self.status.didPass()

    def didDiff(self):
        """
        return boolean for a differential tester failure
        see util.TestStatus for more information
        """
        return self.status.didDiff()

    def isInitialized(self):
        """
        return boolean for tester in an initialization status
        see util.TestStatus for more information
        """
        return self.status.isInitialized()

    def isPending(self):
        """
        return boolean for tester in a pending status
        see util.TestStatus for more information
        """
        return self.status.isPending()

    def isFinished(self):
        """
        return boolean for tester no longer pending
        see util.TestStatus for more information
        """
        return self.status.isFinished()

    def isSkipped(self):
        """
        return boolean for tester being reported as skipped
        see util.TestStatus for more information
        """
        return self.status.isSkipped()

    def isSilent(self):
        """
        return boolean for tester being skipped and not reported
        see util.TestStatus for more information
        """
        return self.status.isSilent()

    def isDeleted(self):
        """
        return boolean for tester being skipped and not reported due to
        internal deletion status
        see util.TestStatus for more information
        """
        return self.status.isDeleted()

    def isQueued(self):
        """
        return boolean for tester in a queued status
        see util.TestStatus for more information
        """
        return self.status.isQueued()

    def isWaiting(self):
        """
        return boolean for tester awaiting process results
        """
        return self.status.isWaiting()

    def getCheckInput(self):
        return self.check_input

    def setValgrindMode(self, mode):
        """ Increase the alloted time for tests when running with the valgrind option """
        if mode == 'NORMAL':
            self.specs['max_time'] = self.specs['max_time'] * 2
        elif mode == 'HEAVY':
            self.specs['max_time'] = self.specs['max_time'] * 6

    def checkRunnable(self, options):
        """
        Derived method to return tuple if this tester should be executed or not.

        The tuple should be structured as (boolean, 'reason'). If false, and the
        reason is left blank, this tester will be treated as silent (no status
        will be printed and will not be counted among the skipped tests).
        """
        return True

    def shouldExecute(self):
        """
        return boolean for tester allowed to execute its command
        see .getCommand for more information
        """
        return self.should_execute

    def prepare(self, options):
        """
        Method which is called prior to running the test. It can be used to cleanup files
        or do other preparations before the tester is run.
        """
        return

    def getThreads(self, options):
        """ return number of threads to use for this tester """
        return 1

    def getProcs(self, options):
        """ return number of processors to use for this tester """
        return 1

    def getSlots(self, options):
        """ return number of slots to use for this tester """
        return self.getThreads(options) * self.getProcs(options)

    def getCommand(self, options):
        """ return the executable command that will be executed by the tester """
        return ''

    def runCommand(self, cmd, cwd, timer, options):
        """
        Helper method for running external (sub)processes as part of the tester's execution.  This
        uses the tester's getCommand and getTestDir methods to run a subprocess.  The timer must
        be the same timer passed to the run method.  Results from running the subprocess is stored
        in the tester's output and exit_code fields.
        """

        cmd = self.getCommand(options)
        cwd = self.getTestDir()

        self.process = None
        try:
            f = TemporaryFile()
            # On Windows, there is an issue with path translation when the command is passed in
            # as a list.
            if platform.system() == "Windows":
                process = subprocess.Popen(cmd,stdout=f,stderr=f,close_fds=False, shell=True, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP, cwd=cwd)
            else:
                process = subprocess.Popen(cmd,stdout=f,stderr=f,close_fds=False, shell=True, preexec_fn=os.setsid, cwd=cwd)
        except:
            print("Error in launching a new task", cmd)
            raise

        self.process = process
        self.outfile = f

        timer.start()
        process.wait()
        timer.stop()

        self.exit_code = process.poll()

        # store the contents of output, and close the file
        self.std_out = util.readOutput(self.outfile, options)
        self.outfile.close()

    def killCommand(self):
        """
        Kills any currently executing process started by the runCommand method.
        """
        if self.process is not None:
            try:
                if platform.system() == "Windows":
                    self.process.terminate()
                else:
                    pgid = os.getpgid(self.process.pid)
                    os.killpg(pgid, SIGTERM)
            except OSError: # Process already terminated
                pass

    def run(self, timer, options):
        """
        This is a method that is the tester's main execution code.  Subclasses can override this
        method with custom code relevant to their specific testing needs.  By default this method
        calls runCommand.  runCommand is provided as a helper for running (external) subprocesses
        as part of the tester's execution and should be the *only* way subprocesses are executed
        if needed. The run method is responsible to call the start+stop methods on timer to record
        the time taken to run the actual test.  start+stop can be called multiple times.
        """
        cmd = self.getCommand(options)
        cwd = self.getTestDir()

        self.runCommand(cmd, cwd, timer, options)

    def processResultsCommand(self, moose_dir, options):
        """ method to return the commands (list) used for processing results """
        return []

    def processResults(self, moose_dir, options, output):
        """ method to process the results of a finished tester """
        return

    def hasRedirectedOutput(self, options):
        """ return bool on tester having redirected output """
        return (self.specs.isValid('redirect_output') and self.specs['redirect_output'] == True and self.getProcs(options) > 1)

    def getRedirectedOutputFiles(self, options):
        """ return a list of redirected output """
        return [os.path.join(self.getTestDir(), self.name() + '.processor.{}'.format(p)) for p in xrange(self.getProcs(options))]

    def addCaveats(self, *kwargs):
        """ Add caveat(s) which will be displayed with the final test status """
        self.__caveats.update(kwargs)
        return self.getCaveats()

    def getCaveats(self):
        """ Return caveats accumalted by this tester """
        return self.__caveats

    def checkRunnableBase(self, options):
        """
        Method to check for caveats that would prevent this tester from
        executing correctly (or not at all).

        DO NOT override this method. Instead, see .checkRunnable()
        """
        reasons = {}
        checks = options._checks

        tag_match = False
        for t in self.tags:
            if t in options.runtags:
                tag_match = True
                break
        if len(options.runtags) > 0 and not tag_match:
            self.setStatus('no tag', self.bucket_silent)
            return False

        # If something has already deemed this test a failure or is silent, return now
        if self.didFail() or self.isSilent():
            return False

        # If --dry-run set the test status to pass and DO NOT return.
        # This will allow additional checks to perform and report tests
        # that would normally be skipped (and return as False).
        if options.dry_run:
            self.success_message = 'DRY RUN'
            self.setStatus(self.success_message, self.bucket_success)

        # Check if we only want to run failed tests
        if options.failed_tests:
            if self.specs['test_name'] not in options._test_list:
                self.setStatus('not failed', self.bucket_silent)
                return False

        # Check if we only want to run syntax tests
        if options.check_input and not self.specs['check_input']:
            self.setStatus('not check_input', self.bucket_silent)
            return False

        # Check if we want to exclude syntax tests
        if options.no_check_input and self.specs['check_input']:
            self.setStatus('is check_input', self.bucket_silent)
            return False

        # Are we running only tests in a specific group?
        if options.group <> 'ALL' and options.group not in self.specs['group']:
            self.setStatus('unmatched group', self.bucket_silent)
            return False
        if options.not_group <> '' and options.not_group in self.specs['group']:
            self.setStatus('unmatched group', self.bucket_silent)
            return False

        # Store regexp for matching tests if --re is used
        if options.reg_exp:
            match_regexp = re.compile(options.reg_exp)

        # If --re then only test matching regexp. Needs to run before other SKIP methods
        # This also needs to be in its own bucket group. We normally print skipped messages.
        # But we do not want to print tests that didn't match regex.
        if options.reg_exp and not match_regexp.search(self.specs['test_name']):
            self.setStatus('silent', self.bucket_silent)
            return False

        # Short circuit method and run this test if we are ignoring all caveats
        if options.ignored_caveats == 'all':
            # Still, we should abide by the derived classes
            return self.checkRunnable(options)

        # Check for deleted tests
        if self.specs.isValid('deleted'):
            reasons['deleted'] = 'deleted ({})'.format(self.specs['deleted'])

        # Skipped by external means (example: TestHarness part2 with --check-input)
        if self.isSkipped():
            reasons['skip'] = self.getStatusMessage()
        # Test is skipped
        elif self.specs.type('skip') is bool and self.specs['skip']:
            # Backwards compatible (no reason)
            reasons['skip'] = 'no reason'
        elif self.specs.type('skip') is not bool and self.specs.isValid('skip'):
            reasons['skip'] = self.specs['skip']
        # If were testing for SCALE_REFINE, then only run tests with a SCALE_REFINE set
        elif (options.scaling) and self.specs['scale_refine'] == 0:
            self.setStatus('silent', self.bucket_silent)
            return False
        # If we're testing with valgrind, then skip tests that require parallel or threads or don't meet the valgrind setting
        elif options.valgrind_mode != '':
            tmp_reason = ''
            if self.specs['valgrind'].upper() == 'NONE':
                tmp_reason = 'Valgrind==NONE'
            elif self.specs['valgrind'].upper() == 'HEAVY' and options.valgrind_mode.upper() == 'NORMAL':
                tmp_reason = 'Valgrind==HEAVY'
            elif self.specs['min_parallel'] > 1 or self.specs['min_threads'] > 1:
                tmp_reason = 'Valgrind requires serial'
            if tmp_reason != '':
                reasons['valgrind'] = tmp_reason
        # If we're running in recover mode skip tests that have recover = false
        elif options.enable_recover and self.specs['recover'] == False:
            reasons['recover'] = 'NO RECOVER'

        # Check for PETSc versions
        (petsc_status, logic_reason, petsc_version) = util.checkPetscVersion(checks, self.specs)
        if not petsc_status:
            reasons['petsc_version'] = 'using PETSc ' + str(checks['petsc_version']) + ' REQ: ' + logic_reason + ' ' + petsc_version

        # Check for SLEPc versions
        (slepc_status, logic_reason, slepc_version) = util.checkSlepcVersion(checks, self.specs)
        if not slepc_status and len(self.specs['slepc_version']) != 0:
            if slepc_version != None:
                reasons['slepc_version'] = 'using SLEPc ' + str(checks['slepc_version']) + ' REQ: ' + logic_reason + ' ' + slepc_version
            elif slepc_version == None:
                reasons['slepc_version'] = 'SLEPc is not installed'

        # PETSc and SLEPc is being explicitly checked above
        local_checks = ['platform', 'compiler', 'mesh_mode', 'method', 'library_mode', 'dtk', 'unique_ids', 'vtk', 'tecplot', \
                        'petsc_debug', 'curl', 'tbb', 'superlu', 'cxx11', 'asio', 'unique_id', 'slepc', 'petsc_version_release', 'boost']
        for check in local_checks:
            test_platforms = set()
            operator_display = '!='
            inverse_set = False
            for x in self.specs[check]:
                if x[0] == '!':
                    if inverse_set:
                        reasons[check] = 'Multiple Negation Unsupported'
                    inverse_set = True
                    operator_display = '=='
                    x = x[1:] # Strip off the !
                x_upper = x.upper()
                if x_upper in test_platforms:
                    reasons[x_upper] = 'Duplicate Entry or Negative of Existing Entry'
                test_platforms.add(x.upper())

            match_found = len(test_platforms.intersection(checks[check])) > 0
            # Either we didn't find the match when we were using normal "include" logic
            # or we did find the match when we wanted to exclude it
            if inverse_set == match_found:
                reasons[check] = re.sub(r'\[|\]', '', check).upper() + operator_display + ', '.join(test_platforms)

        # Check for heavy tests
        if options.all_tests or options.heavy_tests:
            if not self.specs['heavy'] and options.heavy_tests:
                reasons['heavy'] = 'NOT HEAVY'
        elif self.specs['heavy']:
            reasons['heavy'] = 'HEAVY'

        # There should only be one entry in self.specs['dof_id_bytes']
        for x in self.specs['dof_id_bytes']:
            if x != 'ALL' and not x in checks['dof_id_bytes']:
                reasons['dof_id_bytes'] = '--with-dof-id-bytes!=' + x

        # Check to make sure depend files exist
        for file in self.specs['depend_files']:
            if not os.path.isfile(os.path.join(self.specs['base_dir'], file)):
                reasons['depend_files'] = 'DEPEND FILES'

        # We calculate the exe_objects only if we need them
        if self.specs["required_objects"] and checks["exe_objects"] is None:
            checks["exe_objects"] = util.getExeObjects(self.specs["executable"])

        # Check to see if we have the required object names
        for var in self.specs['required_objects']:
            if var not in checks["exe_objects"]:
                reasons['required_objects'] = '%s not found in executable' % var
                break

        # Check to make sure required submodules are initialized
        for var in self.specs['required_submodule']:
            if var not in checks["submodules"]:
                reasons['required_submodule'] = '%s submodule not initialized' % var

        # Check to make sure environment variable exists
        for var in self.specs['env_vars']:
            if not os.environ.has_key(var):
                reasons['env_vars'] = 'ENV VAR NOT SET'

        # Check for display
        if self.specs['display_required'] and not os.getenv('DISPLAY', False):
            reasons['display_required'] = 'NO DISPLAY'

        # Check for queueing
        if (not self.specs['queue_scheduler'] or not self.shouldExecute()) \
           and options.queueing:
            reasons['queue_scheduler'] = 'queue not supported'

        # Remove any matching user supplied caveats from accumulated checkRunnable caveats that
        # would normally produce a skipped test.
        caveat_list = set()
        if options.ignored_caveats:
            caveat_list = set([x.lower() for x in options.ignored_caveats.split()])

        if len(set(reasons.keys()) - caveat_list) > 0:
            tmp_reason = []
            for key, value in reasons.iteritems():
                if key.lower() not in caveat_list:
                    tmp_reason.append(value)

            # Format joined reason to better fit on the screen
            if len(', '.join(tmp_reason)) >= util.TERM_COLS - (len(self.specs['test_name'])+21):
                flat_reason = (', '.join(tmp_reason))[:(util.TERM_COLS - (len(self.specs['test_name'])+24))] + '...'
            else:
                flat_reason = ', '.join(tmp_reason)

            # If the test is deleted we still need to treat this differently
            if 'deleted' in reasons.keys():
                self.setStatus(flat_reason, self.bucket_deleted)
            else:
                self.setStatus(flat_reason, self.bucket_skip)
            return False

        # Check the return values of the derived classes
        self._runnable = self.checkRunnable(options)
        return self._runnable
