#!/usr/bin/env python
import sys
import unittest
import vtk
from PyQt5 import QtWidgets, QtCore

from peacock.ExodusViewer.ExodusViewer import main
from peacock.utils import Testing, qtutils
from mooseutils import message


class TestExodusViewer(Testing.PeacockImageTestCase):
    """
    Testing for ExodusViewer.

    TODO: There is a rendering artifact that shows up in these tests on linux,
          so the imagediffs are not performed.
    """

    #: QApplication: The main App for QT, this must be static to work correctly.
    qapp = QtWidgets.QApplication(sys.argv)

    #: str: The filename to load.
    _filename = Testing.get_chigger_input('mug_blocks_out.e')

    def setUp(self):
        """
        Loads an Exodus file in the VTKWindowWidget object using a structure similar to the ExodusViewer widget.
        """
        message.MOOSE_TESTING_MODE = True
        qtutils.setAppInformation("peacock_exodusviewer")

        settings = QtCore.QSettings()
        settings.clear()
        settings.sync()

        self._widget = main(size=[400,400])
        self._widget.onSetFilenames([self._filename])

        # Start with 'diffused' variable
        self._widget.currentWidget().VariablePlugin.VariableList.setCurrentIndex(2)
        self._widget.currentWidget().VariablePlugin.VariableList.currentIndexChanged.emit(2)


    def write(self, filename):
        """
        Overload the write method.
        """
        self._widget.currentWidget().OutputPlugin.write.emit(filename)

    def testInitial(self):
        """
        Test initial.
        """
        if sys.platform == 'darwin':
            self.assertImage('testInitial.png')
        self.assertFalse(self._widget.cornerWidget().CloseButton.isEnabled())
        self.assertEqual(self._widget.tabText(self._widget.currentIndex()), 'Results')

    def testCloneClose(self):
        """
        Test clone button works.
        """
        self._widget.cornerWidget().clone.emit()
        self._widget.currentWidget().VariablePlugin.VariableList.setCurrentIndex(2)
        self._widget.currentWidget().VariablePlugin.VariableList.currentIndexChanged.emit(2)
        self.assertEqual(self._widget.count(), 2)
        self.assertEqual(self._widget.tabText(self._widget.currentIndex()), 'Results (2)')
        self.assertTrue(self._widget.cornerWidget().CloseButton.isEnabled())

        if sys.platform == 'darwin':
            self.assertImage('testInitial.png')

        # Change camera on cloned tab
        camera = vtk.vtkCamera()
        camera.SetViewUp(-0.7786, 0.2277, 0.5847)
        camera.SetPosition(9.2960, -0.4218, 12.6685)
        camera.SetFocalPoint(0.0000, 0.0000, 0.1250)
        self._widget.currentWidget().VTKWindowPlugin.onCameraChanged(camera)
        if sys.platform == 'darwin':
            self.assertImage('testClone.png')

        # Switch to first tab
        self._widget.setCurrentIndex(0)
        self.assertEqual(self._widget.tabText(self._widget.currentIndex()), 'Results')
        if sys.platform == 'darwin':
            self.assertImage('testInitial.png')

        # Close the first tab
        self._widget.cornerWidget().close.emit()
        self.assertEqual(self._widget.count(), 1)
        self.assertEqual(self._widget.tabText(self._widget.currentIndex()), 'Results (2)')
        self.assertFalse(self._widget.cornerWidget().CloseButton.isEnabled())

    def testMeshState(self):
        """
        Test that the state of the mesh widget after changing files.
        """
        f0 = Testing.get_chigger_input('diffusion_1.e')
        f1 = Testing.get_chigger_input('diffusion_2.e')
        self._widget.onSetFilenames([f0, f1])
        mesh = self._widget.currentWidget().MeshPlugin
        fp = self._widget.currentWidget().FilePlugin
        fp._callbackAvailableFiles(0)
        mesh.ViewMeshToggle.setChecked(False)
        mesh.ScaleX.setValue(.9)
        mesh.ScaleY.setValue(.8)
        mesh.ScaleZ.setValue(.7)
        mesh.Representation.setCurrentIndex(1)
        mesh.DisplacementToggle.setChecked(True)
        mesh.DisplacementMagnitude.setValue(2.0)
        self.assertImage('testDiffusion1.png')

        fp._callbackAvailableFiles(1)
        # had a case where switching files that had the same variable name
        # disabled the entire widget. Couldn't reproduce it with just the MeshPlugin
        # unit tests.
        self.assertEqual(mesh.isEnabled(), True)
        self.assertEqual(mesh.ViewMeshToggle.isEnabled(), True)

        mesh.ViewMeshToggle.setChecked(True)
        mesh.ScaleX.setValue(.7)
        mesh.ScaleY.setValue(.9)
        mesh.ScaleZ.setValue(.8)
        mesh.DisplacementToggle.setChecked(False)
        mesh.DisplacementMagnitude.setValue(1.5)
        self.assertImage('testDiffusion2.png')

        fp._callbackAvailableFiles(0)
        self.assertEqual(mesh.isEnabled(), True)
        self.assertEqual(mesh.ViewMeshToggle.isEnabled(), False) # not enabled for wireframe
        self.assertEqual(mesh.ViewMeshToggle.isChecked(), False)
        self.assertEqual(mesh.ScaleX.value(), .9)
        self.assertEqual(mesh.ScaleY.value(), .8)
        self.assertEqual(mesh.ScaleZ.value(), .7)
        self.assertEqual(mesh.DisplacementToggle.isChecked(), True)
        self.assertEqual(mesh.DisplacementMagnitude.value(), 2.0)
        self.assertImage('testDiffusion1.png')

        fp._callbackAvailableFiles(1)
        self.assertEqual(mesh.isEnabled(), True)
        self.assertEqual(mesh.ViewMeshToggle.isEnabled(), True)
        self.assertEqual(mesh.ViewMeshToggle.isChecked(), True)
        self.assertEqual(mesh.ScaleX.value(), .7)
        self.assertEqual(mesh.ScaleY.value(), .9)
        self.assertEqual(mesh.ScaleZ.value(), .8)
        self.assertEqual(mesh.DisplacementToggle.isChecked(), False)
        self.assertEqual(mesh.DisplacementMagnitude.value(), 1.5)
        self.assertImage('testDiffusion2.png')

    def testPrefs(self):

        settings = QtCore.QSettings()
        settings.setValue("exodus/defaultColorMap", "magma")
        settings.sync()
        self._widget.cornerWidget().clone.emit()
        self.assertEqual(self._widget.preferencesWidget().count(), 6)
        self.assertEqual(self._widget.currentWidget().VariablePlugin.ColorMapList.currentText(), "magma")

        settings.setValue("exodus/defaultColorMap", "default")
        settings.sync()

        self._widget.cornerWidget().clone.emit()
        self.assertEqual(self._widget.currentWidget().VariablePlugin.ColorMapList.currentText(), "default")


if __name__ == '__main__':
    unittest.main(module=__name__, verbosity=2)
