#!/usr/bin/env python
#pylint: disable=missing-docstring
####################################################################################################
#                                    DO NOT MODIFY THIS HEADER                                     #
#                   MOOSE - Multiphysics Object Oriented Simulation Environment                    #
#                                                                                                  #
#                              (c) 2010 Battelle Energy Alliance, LLC                              #
#                                       ALL RIGHTS RESERVED                                        #
#                                                                                                  #
#                            Prepared by Battelle Energy Alliance, LLC                             #
#                               Under Contract No. DE-AC07-05ID14517                               #
#                               With the U. S. Department of Energy                                #
#                                                                                                  #
#                               See COPYRIGHT for full restrictions                                #
####################################################################################################

import os
import unittest
import subprocess
import MooseDocs

class TestPresentation(unittest.TestCase):
    """
    Test that reveal.js slideshow generation is working.
    """
    html_file = os.path.join('examples', 'presentation', 'index.html')
    working_dir = os.getcwd()

    def clean(self):
        """
        Remove the 'presentation.html' file that is used for testing.
        """
        if os.path.exists(self.html_file):
            os.remove(self.html_file)

    def setUp(self):
        """
        Runs prior to each test.
        """
        os.chdir(os.path.join(MooseDocs.ROOT_DIR, 'docs'))
        self.clean()

    def tearDown(self):
        """
        Runs after each test.
        """
        self.clean()
        os.chdir(self.working_dir)

    def testHTML(self):
        """
        Test that slides are generated
        """
        subprocess.check_output(['./moosedocs.py', 'presentation', 'examples/presentation.md'])
        self.assertTrue(os.path.exists(self.html_file))


if __name__ == '__main__':
    unittest.main(verbosity=2)
