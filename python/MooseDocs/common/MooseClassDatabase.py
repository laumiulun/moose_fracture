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
#pylint: enable=missing-docstring
import os
import re
import collections
import MooseDocs

class MooseClassDatabase(object):
    """
    Search the entire repository for *.h files and a look for class definitions.

    This stores the information in a dict() with class name as the key. The value is a tuple
    of ClassInfo named tuples that contain the full path to the file and the remote location
    given the url supplied in the constructor.
    """
    DEFINITION_RE = re.compile(r'class\s*(?P<class>\w+)\b[^;]')
    ClassInfo = collections.namedtuple('ClassInfo', 'filename remote')

    def __init__(self, repo_url):

        self.__repo_url = repo_url.rstrip('/')
        self.__definitions = dict()

        exclude = [os.path.join('moose', 'libmesh'), 'libmesh', '.git', '.lib']
        for base, _, files in os.walk(MooseDocs.ROOT_DIR, topdown=False):
            if any([base.startswith(os.path.join(MooseDocs.ROOT_DIR, sub)) for sub in exclude]):
                continue

            for fname in files:
                if fname.endswith('.h'):
                    full_file = os.path.join(base, fname)
                    self.__search(full_file)

    def __getitem__(self, value):
        """
        Return the dict of class and associated filenames.
        """
        if value in self.__definitions:
            return self.__definitions[value]

    def __contains__(self, value):
        """
        Return if the key exists.
        """
        return value in self.__definitions

    def __search(self, header):
        """
        Search header file for C++ class definition.
        """
        with open(header, 'r') as fid:
            content = fid.read()

        for match in re.finditer(self.DEFINITION_RE, content):
            hfile = header.replace(MooseDocs.ROOT_DIR, '')
            remote = '{}/{}'.format(self.__repo_url, hfile.lstrip('/'))
            store = [self.ClassInfo(filename=hfile, remote=remote)]

            src = header.replace('/include/', '/src/')[:-2] + '.C'
            if os.path.exists(src):
                cfile = src.replace(MooseDocs.ROOT_DIR, '')
                remote = '{}/{}'.format(self.__repo_url, cfile.lstrip('/'))
                store.append(self.ClassInfo(filename=cfile, remote=remote))

            self.__definitions[match.group('class')] = tuple(store)
