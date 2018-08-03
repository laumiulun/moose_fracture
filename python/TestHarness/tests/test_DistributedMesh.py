from TestHarnessTestCase import TestHarnessTestCase

class TestHarnessTester(TestHarnessTestCase):
    def testSyntax(self):
        """
        Test for correct operation with distributed mesh tests
        """

        # Verify the distributed mesh test is skipped
        output = self.runExceptionTests('-i', 'mesh_mode_distributed')
        self.assertIn('skipped (MESH_MODE!=DISTRIBUTED)', output)

        # Verify the distributed mesh test is passing when providing --distributed
        # To be acurate, test for OK rather than asserting if 'distributed' is
        # missing from the output.
        output = self.runTests('--distributed', '-i', 'mesh_mode_distributed')
        self.assertRegexpMatches(output, 'test_harness.distributed_mesh.*?OK')
