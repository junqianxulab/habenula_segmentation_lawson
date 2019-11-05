#!/usr/bin/env python

# DEPENDENCE: PySide, nibabel, numpy


import numpy as np
import nibabel

#from PySide.QtCore import (Qt, QRect, QRectF, QSize, QPointF)
#from PySide.QtCore import Qt
#from PySide.QtGui import (QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPainterPath,
#                              QTransform, QPushButton, QTextEdit, QVBoxLayout, QWidget,
#                              QScrollArea, QSizePolicy, QGraphicsScene, QGraphicsView, QSpinBox,
#                                QPen, QColor, QPainter)

from PySide import QtCore
from PySide.QtGui import (QGridLayout, QHBoxLayout, QLabel, QPainterPath,
                          QTransform, QPushButton, QVBoxLayout, QWidget,
                          QGraphicsScene, QGraphicsView, QSpinBox,
                          QPen, QColor, QPainter)
from PySide.QtGui import (QImage, QPixmap)
from PySide.QtGui import (QCursor)
'''
from qt import QtCore
from qt.QtGui import (QGridLayout, QHBoxLayout, QLabel, QPainterPath,
                        QTransform, QPushButton, QVBoxLayout, QWidget,
                        QGraphicsScene, QGraphicsView, QSpinBox,
                        QPen, QColor, QPainter)
from qt.QtGui import (QImage, QPixmap)
from qt.QtGui import (QCursor)
'''

class HbRegionDraw:
    def __init__(self):
        self.points = [(0, 0)] * 4
        self.insert = 0
        self.area = 0.0
        self.cp = (0, 0)

    def push(self, point, position=None):
        increased = True
        if position is not None:
            self.insert = position
        if self.insert >= 4:
            self.insert = 3
            increased = True
        self.points[self.insert] = point
        self.insert += 1
        return increased

    def pop(self):
        if self.insert < 1:
            return False
        self.insert -= 1
        return True

    def calculate_area(self, name):
        a, b, c, d = [ (float(value[0]), float(value[1])) for value in self.points[:4]]
        #TODO: more than 1 d?
        #self.cp = (a[0], (c[1]-b[1])/(c[0]-b[0])*(a[0]-b[0]) + b[1])
        #self.cp = ( (c[0]-b[0])/(c[1]-b[1])*(a[1]-b[1]) + b[0], a[1] )
        self.cp = ( (c[0]-b[0])/(c[1]-b[1])*(a[1]-c[1]) + c[0], a[1] )
        self.area = 0.0

        if False:
            # Two line segment
            ab = np.sqrt( (b[0]-a[0])**2 + (b[1]-a[1])**2 )
            cpd = np.sqrt( (d[0]-self.cp[0])**2 + (d[1]-self.cp[1])**2 )
            self.area = (ab*cpd)/2

        cp = self.cp
        vA = (a[0]-cp[0], a[1]-cp[1])
        vB = (b[0]-cp[0], b[1]-cp[1])
        vD = (d[0]-cp[0], d[1]-cp[1])

        if name == 'HbRight':
            sign = 1
        else:
            sign = -1

        if vD[0]*vD[0]+vD[1]*vD[1] >= vA[0]*vA[0]+vA[1]*vA[1]:
            rotD = (np.sqrt(vD[0]*vD[0]+vD[1]*vD[1]), 0)
            coefA = rotD[0]
            angleADlad = sign*angleBw2Vectors(vA, vD)
            rotA = rotateVector(vA, angleADlad)
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotA[0]*rotA[0])*rotA[1]*rotA[1])
            eccentricAnomalyADlad = np.arctan(coefA/coefB*np.tan(angleADlad))
            self.area += (np.abs(eccentricAnomalyADlad)*coefA*coefB/2.0)
        else:
            rotA = (np.sqrt(vA[0]*vA[0]+vA[1]*vA[1]), 0)
            coefA = rotA[0]
            angleDAlad = sign*angleBw2Vectors(vD, vA)
            rotD = vD
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotD[0]*rotD[0])*rotD[1]*rotD[1])
            eccentricAnomalyDAlad = np.arctan(coefA/coefB*np.tan(angleDAlad))
            self.area += (np.abs(eccentricAnomalyDAlad)*coefA*coefB/2.0)

        if vD[0]*vD[0]+vD[1]*vD[1] >= vB[0]*vB[0]+vB[1]*vB[1]:
            rotD = (np.sqrt(vD[0]*vD[0]+vD[1]*vD[1]), 0)
            coefA = rotD[0]
            angleADlad = sign*angleBw2Vectors(vA, vD)
            angleBDlad = sign*angleBw2Vectors(vD, vB)
            rotB = rotateVector(vB, angleADlad)
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotB[0]*rotB[0])*rotB[1]*rotB[1])
            eccentricAnomalyDBlad = np.arctan(coefA/coefB*np.tan(angleBDlad))
            self.area += (np.abs(eccentricAnomalyDBlad)*coefA*coefB/2.0)
        else:
            rotB = (np.sqrt(vB[0]*vB[0]+vB[1]*vB[1]), 0)
            coefA = rotB[0]
            angleDBlad = sign*angleBw2Vectors(vD, vB)
            angleABlad = sign*angleBw2Vectors(vA, vB)
            rotD = rotateVector(vD, angleABlad)
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotD[0]*rotD[0])*rotD[1]*rotD[1])
            eccentricAnomalyDBlad = np.arctan(coefA/coefB*np.tan(angleDBlad))
            self.area += (np.abs(eccentricAnomalyDBlad)*coefA*coefB/2.0)

        return self.area

    def get_points(self):
        return self.points

    def write(self, fout):
        for point in self.points:
            fout.write('%s,%s,' % (point[0], point[1]))
        fout.write('%s,%s,%s,%s\n' % (self.cp[0], self.cp[1], self.insert, self.area))

    def read(self, words):
        fwords = [float(word) for word in words]
        self.points = [(fwords[i], fwords[i+1]) for i in range(0, 8, 2)]
        self.cp = (fwords[8], fwords[9])
        self.insert = int(fwords[10])
        self.area = fwords[11]

def writeHb(filename, hbRight, hbLeft):
    fout = open(filename, 'w')
    fout.write('Right\n')
    for key in hbRight:
        if type(key) == type(''):
            continue
        fout.write('%s,' % key)
        hbRight[key].write(fout)
    fout.write('Left\n')
    for key in hbLeft:
        if type(key) == type(''):
            continue
        fout.write('%s,' % key)
        hbLeft[key].write(fout)
    fout.write('End\n')
    fout.close()

def readHb(filename, hbRight, hbLeft, keytype='int'):
    try:
        fin = open(filename)
    except:
        print 'Failed to read %s file' % filename
        return

    for line in fin.readlines():
        if line == 'Right\n':
            hbDict = hbRight
            continue
        if line == 'Left\n':
            hbDict = hbLeft
            continue
        if line == 'End\n':
            break

        words = line[:-1].split(',')
        key = eval('%s(%s)' % (keytype, words[0]))
        value = HbRegionDraw()
        value.read(words[1:])
        hbDict[key] = value

def resetIntensity(image, iMin='auto', iMax='auto', iRange=255):
    if iMin == 'auto':
        iMin = image.min()
    if iMax == 'auto':
        iMax = image.max()

    scaledImage = np.array(image, dtype=float)
    scaledImage[ image < iMin ] = iMin
    scaledImage[ image > iMax ] = iMax

    scaledImage -= iMin
    scaledImage *= iRange
    scaledImage /= (iMax-iMin)
    scaledImage[ scaledImage > iRange ] = iRange
    return scaledImage

def angleBw2Vectors(vector1, vector2):
    divider1 = np.sqrt(vector1[0]**2+vector1[1]**2)
    divider2 = np.sqrt(vector2[0]**2+vector2[1]**2)
    return np.arccos( vector1[0]/divider1*vector2[0]/divider2 + vector1[1]/divider1*vector2[1]/divider2 )

def rotateVector(vector, angle):
    return (vector[0]*np.cos(angle)-vector[1]*np.sin(angle), vector[0]*np.sin(angle)+vector[1]*np.cos(angle))

def convertImageRGB(image):
    shape = list(image.shape)
    if len(shape) != 3:
        print 'image should be 3 dimensional.'
        return False

    shape.append(3)

    if (image > 255).any():
        scaledImage = resetIntensity(image)
    else:
        scaledImage = image
    imageRGB = np.ndarray(shape, dtype=np.int8)
    imageRGB[:,:,:,0] = imageRGB[:,:,:,1] = imageRGB[:,:,:,2] = scaledImage[:]
    return imageRGB

class LoadImage:
    def __init__(self):
        self.zooms = (1,1,1)

    def readNifti1(self, filename):
        img = nibabel.load(filename)
        hdr = img.get_header()
        self.zooms = hdr.get_zooms()
        self.data = img.get_data()

class MyGraphicsScene(QGraphicsScene):
    def __init__(self, imageViewer, viewIndex, parent=None):
        super(MyGraphicsScene, self).__init__(parent)

        self.imageViewer = imageViewer
        self.viewIndex = viewIndex
        self.zoomFactors = imageViewer.zoomFactors

        self.currentPos = None
        self.initialPos = None
        self.pressMiddle = False
        self.pressLeft = False
        self.pressRight = False

    def keyPressEvent(self, event):
        viewIndex = 1
        key = event.key()
        if key == QtCore.Qt.Key_Up:
            self.imageViewer.moveIJK(viewIndex, 1)
        elif key == QtCore.Qt.Key_Down:
            self.imageViewer.moveIJK(viewIndex, -1)
        elif key == QtCore.Qt.Key_Minus:
            self.imageViewer.magZoomFactors(-0.2)
        elif key == QtCore.Qt.Key_Plus:
            self.imageViewer.magZoomFactors(+0.2)
        elif key == QtCore.Qt.Key_Equal:
            self.imageViewer.magZoomFactors(+0.2)

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.MiddleButton:
            self.currentPos = event.scenePos()
            self.pressMiddle = True

        elif event.button() == QtCore.Qt.LeftButton:
            self.currentPos = event.scenePos()
            self.initialPos = event.scenePos()
            self.pressLeft = True

        elif event.button() == QtCore.Qt.RightButton:
            x, y = event.scenePos().toTuple()
            self.imageViewer.mousePressEventRight(self.viewIndex, x, y)

    def mouseMoveEvent(self, event):
        if self.pressMiddle or self.pressLeft:
            currentPos = event.scenePos()
            movingPos = self.currentPos - currentPos
            rect = self.sceneRect()
            topLeft = rect.topLeft()
            rect.moveTo(topLeft + movingPos)
            self.setSceneRect(rect)
            self.currentPos = currentPos + movingPos

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.MiddleButton:
            self.pressMiddle = False
            self.currentPos = None

        elif event.button() == QtCore.Qt.LeftButton:
            if self.initialPos == event.scenePos():
                x, y = event.scenePos().toTuple()
                self.imageViewer.mousePressEventLeft(self.viewIndex, x, y)
            self.pressMiddle = False
            self.currentPos = None
            self.initialPos = None

    def savePng(self, filename, filetype='PNG'):
        rect = self.sceneRect()
        pix = QPixmap(rect.width(), rect.height())
        painter = QPainter(pix)
        self.render(painter, painter.window(), rect)
        pix.save(filename, filetype)
        del painter
        del pix

class ImageViewer(QWidget):
    def __init__(self, filename=None, name=None, parent=None):
        super(ImageViewer, self).__init__(parent)
        self.imageRaw = np.zeros((200, 200, 200))
        self.imageRGB = np.zeros((200, 200, 200, 3))
        self.viewRange = [[100, 100], [100, 100], [100, 100]]
        self.zoomFactors = [3.0, 3.0, 3.0]
        self.iMin = 0
        self.iMax = 1000
        self.ijk = [0,0,0]

        self.name = 'untitled'
        self.initHbRL()
        self.mouseCursorCross = False

        ## Upper buttons
        self.upperRow = QHBoxLayout()
        self.drawRightHb = QPushButton('RightHb', self)
        self.drawLeftHb = QPushButton('LeftHb', self)
        self.drawRightHb.setCheckable(True)
        self.drawLeftHb.setCheckable(True)
        self.upperRow.addWidget(self.drawRightHb)
        self.upperRow.addWidget(self.drawLeftHb)
        QtCore.QObject.connect(self.drawRightHb, QtCore.SIGNAL('clicked()'), self.clickRightHb)
        QtCore.QObject.connect(self.drawLeftHb, QtCore.SIGNAL('clicked()'), self.clickLeftHb)

        ## View Layout
        self.scenes = [MyGraphicsScene(self, 0), MyGraphicsScene(self, 1), MyGraphicsScene(self, 2)]
        self.views = [QGraphicsView(self.scenes[0]), QGraphicsView(self.scenes[1]), QGraphicsView(self.scenes[2])]
        self.sceneItems = [None, None, None]

        #Scene 4
        view4 = QGridLayout()
        labelX = QLabel('X: ')
        labelY = QLabel('Y: ')
        labelZ = QLabel('Z: ')
        self.spin = [QSpinBox(), QSpinBox(), QSpinBox()]
        labelMinIntensity = QLabel('Min: ')
        labelMaxIntensity = QLabel('Max: ')
        labelZoom = QLabel('Zoom: ')
        self.spinMin = QSpinBox()
        self.spinMax = QSpinBox()
        self.spinZoom = QSpinBox()
        self.spinMin.setMaximum(10000)
        self.spinMax.setMaximum(10000)
        self.spinMin.setValue(self.iMin)
        self.spinMax.setValue(self.iMax)
        self.spinZoom.setValue(self.zoomFactors[0])
        self.spinMin.setKeyboardTracking(False)
        self.spinMin.valueChanged.connect(self.setIntensity)
        self.spinMax.setKeyboardTracking(False)
        self.spinMax.valueChanged.connect(self.setIntensityMax)
        self.spinZoom.valueChanged.connect(self.setZoomFactor)
        self.spinZoom.setKeyboardTracking(False)
        labelAreaR = QLabel('Right: ')
        labelAreaL = QLabel('Left: ')
        self.labelAreaValueR = QLabel()
        self.labelAreaValueL = QLabel()

        view4.addWidget(labelX, 0, 0, 2, 1)
        view4.addWidget(labelY, 1, 0, 2, 1)
        view4.addWidget(labelZ, 2, 0, 2, 1)
        view4.addWidget(labelMinIntensity, 4, 0, 2, 1)
        view4.addWidget(labelMaxIntensity, 5, 0, 2, 1)
        view4.addWidget(labelZoom, 6, 0, 2, 1)
        view4.addWidget(labelAreaR, 7, 0)
        view4.addWidget(labelAreaL, 7, 2)
        view4.addWidget(self.spin[0], 0, 2, 2, 1)
        view4.addWidget(self.spin[1], 1, 2, 2, 1)
        view4.addWidget(self.spin[2], 2, 2, 2, 1)
        view4.addWidget(self.spinMin, 4, 2, 2, 1)
        view4.addWidget(self.spinMax, 5, 2, 2, 1)
        view4.addWidget(self.spinZoom, 6, 2, 2, 1)
        view4.addWidget(self.labelAreaValueR, 7, 1)
        view4.addWidget(self.labelAreaValueL, 7, 3)

        self.viewLayout = QGridLayout()
        self.viewLayout.addWidget(self.views[0], 0, 0)
        self.viewLayout.addWidget(self.views[1], 0, 1)
        self.viewLayout.addWidget(self.views[2], 1, 0)
        self.viewLayout.addLayout(view4, 1, 1)

        ## Lower
        self.lowerRow = QHBoxLayout()
        calculateButtonDown = QPushButton('calculateVolume', self)
        saveButtonDown = QPushButton('SaveHb', self)
        loadButtonDown = QPushButton('loadHb', self)
        savepngButtonDown = QPushButton('Save PNG', self)
        self.lowerRow.addWidget(calculateButtonDown)
        self.lowerRow.addWidget(saveButtonDown)
        self.lowerRow.addWidget(loadButtonDown)
        self.lowerRow.addWidget(savepngButtonDown)
        QtCore.QObject.connect(calculateButtonDown, QtCore.SIGNAL('clicked()'), self.getHbAreas)
        QtCore.QObject.connect(saveButtonDown, QtCore.SIGNAL('clicked()'), self.saveHb)
        QtCore.QObject.connect(loadButtonDown, QtCore.SIGNAL('clicked()'), self.loadHb)
        QtCore.QObject.connect(savepngButtonDown, QtCore.SIGNAL('clicked()'), self.saveCoronalSlice)

        ## Layout
        self.mainLayout = QVBoxLayout()
        self.mainLayout.addLayout(self.upperRow)
        self.mainLayout.addLayout(self.viewLayout)
        self.mainLayout.addLayout(self.lowerRow)

        self.setLayout(self.mainLayout)
        self.setWindowTitle("Image Viewer")
        self.resize(600, 700)
        self.show()

        if filename:
            self.loadNifti1(filename, name=name)

    def loadNifti1(self, filename, name=None):
        self.loadImage = LoadImage()
        self.loadImage.readNifti1(filename)
        if name is None:
            name = filename[:]
            if name[-7:] == '.nii.gz':
                name = name[:-7]
            elif name[-4:] == '.nii':
                name = name[:-4]

            for ch in ['/', '.', '~']:
                pos = name.find(ch)
                while pos > -1:
                    name = '%s_%s' % (name[:pos], name[pos+1:])
                    pos = name.find(ch)


        self.name = name
        self.imageRaw = self.loadImage.data

        self.spinMin.setMaximum(int(self.imageRaw.max()*1.2))
        self.spinMax.setMaximum(int(self.imageRaw.max()*1.2))
        #self.iMin = int(max(0, self.imageRaw.min()*1.2))
        self.iMin = 0
        self.iMax = int(self.imageRaw.max()*0.8)

        self.spinMin.setValue(self.iMin)
        self.spinMax.setValue(self.iMax)

        self.imageRGB = convertImageRGB(resetIntensity(self.imageRaw, iMin=self.iMin, iMax=self.iMax))
        self.setIJK()
        self.drawBaseImage()
        self.setSceneIJK()
        self.showSpinValue()

        self.setMouseCursor(True)

    def showSpinValue(self):
        self.spin[0].setMaximum(self.imageRaw.shape[0])
        self.spin[1].setMaximum(self.imageRaw.shape[1])
        self.spin[2].setMaximum(self.imageRaw.shape[2])
        self.spin[0].setKeyboardTracking(False)
        self.spin[1].setKeyboardTracking(False)
        self.spin[2].setKeyboardTracking(False)
        self.spin[0].setValue(self.ijk[0])
        self.spin[1].setValue(self.ijk[1])
        self.spin[2].setValue(self.ijk[2])
        self.spin[0].valueChanged.connect(self.resetI)
        self.spin[1].valueChanged.connect(self.resetJ)
        self.spin[2].valueChanged.connect(self.resetK)

    def resetI(self, i):
        self.ijk[0] = i
        self.drawBaseImage([0])
    def resetJ(self, j):
        self.ijk[1] = j
        self.drawBaseImage([1])
        self.labelAreaValueR.clear()
        self.labelAreaValueL.clear()
    def resetK(self, k):
        self.ijk[2] = k
        self.drawBaseImage([2])

    def initHbRL(self):
        self.hbRight = {'name':'HbRight', 'color':'red', 'images':[], 'points':[]}
        self.hbLeft = {'name':'HbLeft', 'color':'blue', 'images':[], 'points':[]}

    def clickRightHb(self):
        if self.drawLeftHb.isChecked():
            self.drawLeftHb.setChecked(False)

    def clickLeftHb(self):
        if self.drawRightHb.isChecked():
            self.drawRightHb.setChecked(False)

    def saveHb(self, filename=None):
        if filename is None:
            filename = '%s_lawson_hb.csv' % self.name
        writeHb(filename, self.hbRight, self.hbLeft)
        with open(filename, 'a') as fout:
            volumes = self.getHbAreas()
            fout.write('%s,%s\n' % (volumes[0], volumes[1]) )

        png_filename = filename[:-4]

    def saveCoronalSlice(self, filename=None):
        viewIndex = 1
        if filename is None:
            filename = '%s_lawson_hb' % self.name
        filename += str('%03d.png' % self.ijk[viewIndex])
        self.scenes[viewIndex].savePng(filename, 'PNG')

    def loadHb(self, filename=None):
        if filename is None:
            filename = '%s_lawson_hb.csv' % self.name
        self.initHbRL()
        readHb(filename, self.hbRight, self.hbLeft)

    def setMouseCursor(self, cross=False):
        if cross != self.mouseCursorCross:
            self.mouseCursorCross = cross
            if cross:
                #self.setCursor(QCursor(QtCore.Qt.CrossCursor))
                self.setCursor(QtCore.Qt.CrossCursor)
            else:
                #self.setCursor(QCursor(QtCore.Qt.ArrowCursor))
                self.setCursor(QtCore.Qt.ArrowCursor)

    def getZoomFactors(self):
        return self.zoomFactors

    def magZoomFactors(self, ratio):
        self.setZoomFactors([value+ratio for value in self.zoomFactors])
        self.spinZoom.setValue(self.zoomFactors[0])

    def setZoomFactors(self, zoomFactors):
        sceneRects = self.saveScene()
        self.zoomFactors = zoomFactors
        self.drawBaseImage()
        self.restoreScene(sceneRects)

    def setZoomFactor(self, zoomFactor):
        self.setZoomFactors([zoomFactor, zoomFactor, zoomFactor])

    def getIntensity(self):
        return self.iMin, self.iMax

    def setIntensityMax(self, iMax):
        self.setIntensity(iMax=iMax)
    def setIntensity(self, iMin=None, iMax=None):
        if iMin:
            self.iMin = iMin
        if iMax:
            self.iMax = iMax
        self.imageRGB = convertImageRGB(resetIntensity(self.imageRaw, iMin=self.iMin, iMax=self.iMax))
        self.drawBaseImage()

    def mousePressEventLeft(self, viewIndex, x, y):
        if viewIndex == 1:
            if self.drawRightHb.isChecked():
                hbDict = self.hbRight
                labelAreaValue = self.labelAreaValueR
            elif self.drawLeftHb.isChecked():
                hbDict = self.hbLeft
                labelAreaValue = self.labelAreaValueL
            else:
                self.resetIJK(viewIndex, x, y)
                return

            if not hbDict.has_key(self.ijk[viewIndex]):
                hbDict[self.ijk[viewIndex]] = HbRegionDraw()

            if hbDict[self.ijk[viewIndex]].insert == 4:
                for item in hbDict['images']:       # lines
                    self.scenes[viewIndex].removeItem(item)
                hbDict['images'] = []
                item = hbDict['points'].pop(0)       # last point
                self.scenes[viewIndex].removeItem(item)

            i, j = self.point2fijk(viewIndex, x, y)
            if not hbDict[self.ijk[viewIndex]].push((i,j)):
                item = hbDict['points'].pop(0)       # last point
                self.scenes[viewIndex].removeItem(item)
                #self.scenes[viewIndex].removeItem(self.scenes[viewIndex].items()[0])

            x_round, y_round = self.ij2point(viewIndex, i, j)
            hbDict['points'].insert(0, self.scenes[viewIndex].addEllipse(x_round-1, y_round-1, 2, 2, pen=QPen(QColor(hbDict['color']))))

            if hbDict[self.ijk[viewIndex]].insert >= 4:
                area = hbDict[self.ijk[viewIndex]].calculate_area(hbDict['name'])
                print ' %s %s: %s' % (self.ijk[viewIndex], hbDict['name'], area)
                self.drawHb(hbDict)
                #labelAreaValue.setText('%s: %3.2f' % (self.ijk[viewIndex], area))
                labelAreaValue.setText('%3.2f' % (area))
        else:
            self.resetIJK(viewIndex, x, y)

    def mousePressEventRight(self, viewIndex, x, y):
        if viewIndex == 1:
            if self.drawRightHb.isChecked():
                hbDict = self.hbRight
            elif self.drawLeftHb.isChecked():
                hbDict = self.hbLeft
            else:
                return

            if not hbDict.has_key(self.ijk[viewIndex]):
                return

            if hbDict[self.ijk[viewIndex]].insert == 4:
                for item in hbDict['images']:       # lines
                    self.scenes[viewIndex].removeItem(item)
                hbDict['images'] = []
                #for item in self.scenes[viewIndex].items()[:3]:     # lines
                #    self.scenes[viewIndex].removeItem(item)

            if hbDict[self.ijk[viewIndex]].pop():
                item = hbDict['points'].pop(0)       # last point
                self.scenes[viewIndex].removeItem(item)
                #self.scenes[viewIndex].removeItem(self.scenes[viewIndex].items()[0])

    def moveIJK(self, viewIndex, dx):
        self.ijk[viewIndex] += dx
        self.spin[viewIndex].setValue(self.ijk[viewIndex])
        self.drawBaseImage([viewIndex])
        self.labelAreaValueR.clear()
        self.labelAreaValueL.clear()

    def point2fijk(self, viewIndex, x, y):
        indexTable = [[1, 2], [0, 2], [0, 1]]
        xIndex = indexTable[viewIndex][0]
        yIndex = indexTable[viewIndex][1]
        ijkX = float(x)/self.zoomFactors[viewIndex]
        ijkY = self.imageRGB.shape[yIndex] - 1 - float(y)/self.zoomFactors[viewIndex]
        ijkX = np.round(ijkX*2)/2
        ijkY = np.round(ijkY*2)/2
        return (ijkX, ijkY)

    def point2fijk_flt(self, viewIndex, x, y):
        indexTable = [[1, 2], [0, 2], [0, 1]]
        xIndex = indexTable[viewIndex][0]
        yIndex = indexTable[viewIndex][1]
        ijkX = float(x)/self.zoomFactors[viewIndex]
        ijkY = self.imageRGB.shape[yIndex] - 1 - float(y)/self.zoomFactors[viewIndex]
        return (ijkX, ijkY)

    def point2ijk(self, viewIndex, x, y):
        indexTable = [[1, 2], [0, 2], [0, 1]]
        xIndex = indexTable[viewIndex][0]
        yIndex = indexTable[viewIndex][1]
        ijkX = int(x/self.zoomFactors[viewIndex])
        ijkY = self.imageRGB.shape[yIndex] - 1 - int(y/self.zoomFactors[viewIndex])
        return (ijkX, ijkY)

    def ij2point(self, viewIndex, i, j):
        indexTable = [[1, 2], [0, 2], [0, 1]]
        xIndex = indexTable[viewIndex][0]
        yIndex = indexTable[viewIndex][1]
        pointX = i*self.zoomFactors[viewIndex]
        pointY = (self.imageRGB.shape[yIndex] - 1 - j)*self.zoomFactors[viewIndex]
        return (pointX, pointY)

    def resetIJK(self, viewIndex, x, y):
        indexTable = [[1, 2], [0, 2], [0, 1]]
        xIndex = indexTable[viewIndex][0]
        yIndex = indexTable[viewIndex][1]
        self.ijk[xIndex] = int(x/self.zoomFactors[viewIndex])
        self.ijk[yIndex] = self.imageRGB.shape[yIndex] - 1 - int(y/self.zoomFactors[viewIndex])
        self.spin[xIndex].setValue(self.ijk[xIndex])
        self.spin[yIndex].setValue(self.ijk[yIndex])
        self.drawBaseImage(indexTable[viewIndex])
        self.labelAreaValueR.clear()
        self.labelAreaValueL.clear()

    def saveScene(self):
        sceneRects = [0, 0, 0]
        for i in range(3):
            size = self.views[i].size().toTuple()
            x, y, _, _ = self.scenes[i].sceneRect().getRect()
            x += (size[0]-10)/2.0
            y += (size[1]-10)/2.0
            x /= self.zoomFactors[i]
            y /= self.zoomFactors[i]
            sceneRects[i] = [x, y]
        return sceneRects

    def restoreScene(self, sceneRects):
        for i in range(3):
            size = self.views[i].size().toTuple()
            x = sceneRects[i][0]*self.zoomFactors[i]
            y = sceneRects[i][1]*self.zoomFactors[i]
            self.scenes[i].setSceneRect(x-(size[0]-10)/2, y-(size[0]-10)/2, size[0]-10, size[1]-10)

    def setSceneIJK(self, viewIndex=(0,1,2)):
        for i in viewIndex:
            size = self.views[i].size().toTuple()

            xyCenter = [ self.ijk[tmp]*self.zoomFactors[i] for tmp in range(3) if tmp is not i]

            self.scenes[i].setSceneRect(xyCenter[0] - size[0]/2.0, xyCenter[1] - size[1]/2.0, size[0]-10, size[1]-10)

    def setIJK(self, ijk='center'):
        if ijk == 'center':
            self.ijk = [self.imageRGB.shape[0]/2, self.imageRGB.shape[1]/2, self.imageRGB.shape[2]/2]
        else:
            self.ijk = ijk
        self.spin[0].setValue(self.ijk[0])
        self.spin[1].setValue(self.ijk[1])
        self.spin[2].setValue(self.ijk[2])

    def get2D(self, plane):
        if plane == 0:
            return np.rot90(self.imageRGB[self.ijk[0], :, :, :])
        elif plane == 1:
            return np.rot90(self.imageRGB[:, self.ijk[1], :, :])
        elif plane == 2:
            return np.rot90(self.imageRGB[:, :, self.ijk[2], :])

    def drawBaseImage(self, viewIndex=(0,1,2), eraseOthers=True):
        for i in viewIndex:
            imageArray3d = self.get2D(i)
            height, width, bytesPerComponent = imageArray3d.shape
            imageArray = imageArray3d.flatten()
            bytesPerLine = bytesPerComponent * width
            image = QImage(imageArray, width, height, bytesPerLine, QImage.Format_RGB888)
            #image = image.scaled( int(width*self.zoomFactors[i]), int(height*self.zoomFactors[i]), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            image = image.scaled( int(width*self.zoomFactors[i]), int(height*self.zoomFactors[i]), QtCore.Qt.KeepAspectRatio)

            if self.sceneItems[i] is not None:
                #self.scenes[i].removeItem(self.sceneItems[i])
                self.sceneItems[i].setPixmap(QPixmap.fromImage(image))
                if False:
                    items = self.scenes[i].items()
                    for item in items:
                        self.scenes[i].removeItem(item)

                    self.sceneItems[i] = self.scenes[i].addPixmap(QPixmap.fromImage(image))
                    for item in items[1:]:
                        self.scenes[i].additem(item)

                    del items

            else:
                self.sceneItems[i] = self.scenes[i].addPixmap(QPixmap.fromImage(image))

            if eraseOthers:
                for item in self.scenes[i].items():
                    if item != self.sceneItems[i]:
                        self.scenes[i].removeItem(item)

            # TODO: where to put this
            self.drawHbPts(self.hbRight)
            self.drawHb(self.hbRight)
            self.drawHbPts(self.hbLeft)
            self.drawHb(self.hbLeft)

    def drawHbPts(self, hbDict):
        viewIndex = 1
        if not hbDict.has_key(self.ijk[viewIndex]):
            return
        pts = [ self.ij2point(viewIndex, i, j) for (i,j) in hbDict[self.ijk[viewIndex]].get_points() ]
        for i in range(hbDict[self.ijk[viewIndex]].insert):
            hbDict['points'].insert(0, self.scenes[viewIndex].addEllipse(pts[i][0]-1, pts[i][1]-1, 2, 2, pen=QPen(QColor(hbDict['color']))))

    def drawHb(self, hbDict):
        viewIndex = 1
        if not hbDict.has_key(self.ijk[viewIndex]):
            return
        if hbDict[self.ijk[viewIndex]].insert < 4:
            return

        a, b, c, d = [ self.ij2point(viewIndex, i, j) for (i,j) in hbDict[self.ijk[viewIndex]].get_points() ]

        cp = hbDict[self.ijk[viewIndex]].cp
        cp = self.ij2point(viewIndex, cp[0], cp[1])

        vA = (a[0]-cp[0], a[1]-cp[1])
        vB = (b[0]-cp[0], b[1]-cp[1])
        vD = (d[0]-cp[0], d[1]-cp[1])

        arcAD = QPainterPath()
        arcDB = QPainterPath()

        if hbDict['name'] == 'HbRight':
            sign = 1
            startAngle = 0
        else:
            sign = -1
            startAngle = 180

        angleADlad = sign*angleBw2Vectors(vA, vD)
        angleAD = 180/np.pi*angleADlad

        if vD[0]*vD[0]+vD[1]*vD[1] >= vA[0]*vA[0]+vA[1]*vA[1]:
            rotD = (np.sqrt(vD[0]*vD[0]+vD[1]*vD[1]), 0)
            coefA = rotD[0]

            rotA = rotateVector(vA, angleADlad)
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotA[0]*rotA[0])*rotA[1]*rotA[1])
            eccentricAnomalyAD = -180/np.pi*np.arctan(coefA/coefB*np.tan(angleADlad))
            arcAD.moveTo(cp[0]+sign*coefA, cp[1])
            arcAD.arcTo(cp[0]-coefA, cp[1]-coefB, 2*coefA, 2*coefB, startAngle, eccentricAnomalyAD)
        else:
            rotA = (np.sqrt(vA[0]*vA[0]+vA[1]*vA[1]), 0)
            coefA = rotA[0]

            angleDAlad = sign*angleBw2Vectors(vD, vA)
            angleDA = 180/np.pi*angleDAlad

            rotD = vD
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotD[0]*rotD[0])*rotD[1]*rotD[1])
            eccentricAnomalyDA = 180/np.pi*np.arctan(coefA/coefB*np.tan(angleDAlad))
            arcAD.moveTo(cp[0]+sign*coefA, cp[1])
            arcAD.arcTo(cp[0]-coefA, cp[1]-coefB, 2*coefA, 2*coefB, startAngle, eccentricAnomalyDA)
            angleAD = 0.0

        if vD[0]*vD[0]+vD[1]*vD[1] >= vB[0]*vB[0]+vB[1]*vB[1]:
            rotD = (np.sqrt(vD[0]*vD[0]+vD[1]*vD[1]), 0)
            coefA = rotD[0]

            angleBDlad = sign*angleBw2Vectors(vD, vB)
            angleBD = 180/np.pi*angleBDlad

            rotB = rotateVector(vB, angleADlad)
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotB[0]*rotB[0])*rotB[1]*rotB[1])
            eccentricAnomalyDB = 180/np.pi*np.arctan(coefA/coefB*np.tan(angleBDlad))
            arcDB.moveTo(cp[0]+sign*coefA, cp[1])
            arcDB.arcTo(cp[0]-coefA, cp[1]-coefB, 2*coefA, 2*coefB, startAngle, eccentricAnomalyDB)

            angleAB = 180/np.pi*sign*angleBw2Vectors(vA, vD)
        else:
            rotB = (np.sqrt(vB[0]*vB[0]+vB[1]*vB[1]), 0)
            coefA = rotB[0]

            angleABlad = sign*angleBw2Vectors(vA, vB)
            angleAB = 180/np.pi*angleABlad

            angleDBlad = sign*angleBw2Vectors(vD, vB)
            angleDB = 180/np.pi*angleDBlad

            rotD = rotateVector(vD, angleABlad)
            coefB = np.sqrt(coefA*coefA/(coefA*coefA-rotD[0]*rotD[0])*rotD[1]*rotD[1])
            eccentricAnomalyDB = -180/np.pi*np.arctan(coefA/coefB*np.tan(angleDBlad))
            arcDB.moveTo(cp[0]+sign*coefA, cp[1])
            arcDB.arcTo(cp[0]-coefA, cp[1]-coefB, 2*coefA, 2*coefB, startAngle, eccentricAnomalyDB)
            

        imgAC = self.scenes[viewIndex].addLine(a[0], a[1], cp[0], cp[1], pen=QPen(QColor(hbDict['color'])))
        imgBC = self.scenes[viewIndex].addLine(b[0], b[1],  c[0],  c[1], pen=QPen(QColor(hbDict['color'])))
        imgAD = self.scenes[viewIndex].addPath(arcAD, pen=QPen(QColor(hbDict['color'])))
        imgAD.setTransform(QTransform().translate(cp[0], cp[1]).rotate(-angleAD).translate(-cp[0], -cp[1]))
        imgDB = self.scenes[viewIndex].addPath(arcDB, pen=QPen(QColor(hbDict['color'])))
        imgDB.setTransform(QTransform().translate(cp[0], cp[1]).rotate(-angleAB).translate(-cp[0], -cp[1]))
        hbDict['images'] = [imgAD, imgDB, imgBC, imgAC]

    def getHbAreas(self):
        volumes = []
        for hbDict in [self.hbRight, self.hbLeft]:
            keys = hbDict.keys()
            if not keys:
                volumes.append(0.0)
                continue

            print 'Calculate %s Volume...' % hbDict['name']
            unit = self.loadImage.zooms[0] * self.loadImage.zooms[1] * self.loadImage.zooms[2]
            volume = 0.0
            keys.sort()
            for key in keys:
                if type(key) == type(''):
                    continue
                area = hbDict[key].area
                print ' %s: %s (voxel)' % (key, area)
                volume += area
            volume *= unit
            print 'volume = %s (mm^3)' % volume
            volumes.append(volume)
        return volumes

if __name__ == '__main__':
    import sys
    from PySide.QtGui import QApplication

    if len(sys.argv) > 1:
        filename = sys.argv[1]
        if len(sys.argv) > 2:
            name = sys.argv[2]
        else:
            name = None
    else:
        print 'Usage: %s image_filename name_for_output' % sys.argv[0]
        sys.exit()


    app = QApplication(sys.argv)

    imageViewer = ImageViewer(filename=filename, name=name)
    #imageViewer.show()
    #imageViewer.setSceneIJK()

    sys.exit(app.exec_())


