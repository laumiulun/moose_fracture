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
import logging
from markdown.util import etree

LOG = logging.getLogger(__name__)

class MooseCollapsible(object):
    """
    A tool for creating materialize collections.
    """

    def __init__(self):
        self._ul = etree.Element('ul')
        self._ul.set('class', 'collapsible')
        self._ul.set('data-collapsible', 'accordion')
        self.__count = 0

    def __nonzero__(self):
        return self.__count > 0

    def addHeader(self, header):
        """
        Adds a collapsible-header div to the list with the 'moose-group-header' class for showing
        the group divisions of objects and systems.
        """
        self.__count += 1
        div = etree.SubElement(self._ul, 'div')
        div.set('class', 'collapsible-header moose-group-header')
        div.text = header

    def addItem(self, name, desc=None, body=None, id_='moose-collection-li'):
        """
        Add a collapsible item to the list.
        """
        self.__count += 1

        li = etree.SubElement(self._ul, 'li')
        li.set('id', id_)

        header = etree.SubElement(li, 'div')
        header.set('class', 'collapsible-header')

        if desc:
            btn = etree.SubElement(header, 'i')
            btn.set('class', 'material-icons')
            btn.text = 'keyboard_arrow_down'

        name_el = etree.SubElement(header, 'span')
        name_el.set('class', 'moose-collection-name')
        if isinstance(name, (unicode, str)):
            name_el.text = name
        else:
            name_el.append(name)

        if desc:
            desc_el = etree.SubElement(header, 'span')
            desc_el.set('class', 'moose-collection-description')
            desc_el.text = ": " + desc

        if body:
            body_el = etree.SubElement(li, 'div')
            body_el.set('class', 'collapsible-body')
            body_el.text = body

    def element(self):
        """
        Return the unordered list tag.
        """
        return self._ul
