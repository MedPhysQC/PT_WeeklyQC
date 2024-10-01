#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""

Dependencies:

    PYTHON
    - pydicom / dicom
    - numpy
    - datetime
    - matplotlib
    - string


Version control:
    20181128: first version
    20190227: [JK] add deviation from calibration per quadrant
    20190228: [JK] reverse output images upside down, change slc/quadant/metric indexing
    20190306: [JK] fix issue with text on quadrant images
    20220928: [JK] get acquisition time from DICOM tag AcquisitionTime (was: last part of DateTime)


The initial version was created by Stijn van de Schoot at the Amsterdam UMC, location VUmc

"""

__name__ = 'PETQC_lib.py'
__version__ = '20220928'
__author__ = 'Stijn van de Schoot, Maqsood Yaqub'


import matplotlib
matplotlib.use('Agg')


import math
import os
import numpy as np
import datetime
import matplotlib.pyplot as plt
import string

# Import DICOM
try:
    import pydicom as dicom
except ImportError:
    import dicom



class IndexConverter:

    def __init__(self):
        self.maxAlphabet = 26

    def indexToLetter(self, index):
        # Get a letter based on the index
        alphabet = list(string.ascii_uppercase)
        if index < self.maxAlphabet:
            letter = alphabet[index]
        else:
            tmpX = int(index / self.maxAlphabet)
            firstLetter = alphabet[tmpX - 1]
            letter = firstLetter + alphabet[index - (tmpX * self.maxAlphabet)]
        return letter



class ROI:

    def __init__(self, params, dcmFile):
        self.radRoiCM = float(params["ROI_radius"])
        self.threshold = int(params["threshold_value"])
        self.dcmData = dicom.read_file(dcmFile)

        # initiate pixel dictionary
        self.pixelDictionary = {}

        # get pixel spacing from dicom file and translate from mm to cm.
        if 'PixelSpacing' in self.dcmData:
            self.pixelSpacingRowCM = self.dcmData.PixelSpacing[0] / 10
            self.pixelSpacingColumnCM = self.dcmData.PixelSpacing[1] / 10

            if self.pixelSpacingRowCM != self.pixelSpacingColumnCM:
                print('[{0}] :: GetPixelsInROI - Initialisation -- ERROR: The data consists of non-uniform pixels. At this moment, this is not supported by the module!'.format(os.path.basename(__file__)))

        else:
            print('[{0}] :: GetPixelsInROI - Initialisation -- ERROR: The DICOM file does not contain pixel spacing information!'.format(os.path.basename(__file__)))


    def getCentralPixel(self):
        """
        Determine the central pixel of the phantom on the 2D dicom image based on intensity thresholding
        """

        circle_count = 0
        circle_sum_x = 0
        circle_sum_y = 0
        circle_min_x = -1
        circle_max_x = 0

        for y in range(self.dcmData.Rows):
            for x in range(self.dcmData.Columns):
                if self.dcmData.pixel_array[y, x] > self.threshold:
                    circle_count += 1
                    circle_sum_x += x
                    circle_sum_y += y
                    if circle_min_x == -1:
                        circle_min_x = x
                    else:
                        circle_min_x = min(circle_min_x, x)
                    circle_max_x = max(circle_max_x, x)

        if circle_count == 0:
            print('[{0}] :: GetPixelsInROI - getCentralPixel -- ERROR: No phantom detected using the current threshold value!'.format(os.path.basename(__file__)))
            self.Phantom_center_Column = -1
            self.Phantom_center_Row = -1
            self.Phantom_radius = -1
        else:
            self.Phantom_center_Column = round((circle_sum_x / circle_count) * 2) / 2
            self.Phantom_center_Row = round((circle_sum_y / circle_count) * 2) / 2
            self.Phantom_radiusCM = (((circle_max_x - circle_min_x) +1) / 2) * self.pixelSpacingColumnCM

            print('[{0}] :: GetPixelsInROI - getCentralPixel -- The number of pixels detected as being phantom: {1}'.format(os.path.basename(__file__), circle_count))
            print('[{0}] :: GetPixelsInROI - getCentralPixel -- The detected central pixel of the phantom is: x={1}; y={2}'.format(os.path.basename(__file__), self.Phantom_center_Column, self.Phantom_center_Row))
            print('[{0}] :: GetPixelsInROI - getCentralPixel -- The detected radius of the phantom is: {1}'.format(os.path.basename(__file__), self.Phantom_radiusCM))


    def getPixelsInROI(self):
        # initiate/clear pixel dictionary
        self.pixelDictionary = {}
        # translate user-defined ROI to pixels
        radRoiPixels = self.radRoiCM / self.pixelSpacingColumnCM

        # determine bounding box of ROI
        print('[{0}] :: GetPixelsInROI - getPixelsInROI -- Determining the bounding box of the ROI.'.format(os.path.basename(__file__)))
        iLeft = int(math.floor(self.Phantom_center_Column - radRoiPixels))
        iRight = int(math.ceil(self.Phantom_center_Column + radRoiPixels))
        iTop = int(math.floor(self.Phantom_center_Row - radRoiPixels))
        iBottom = int(math.ceil(self.Phantom_center_Row + radRoiPixels))
        print('[{0}] :: GetPixelsInROI - getPixelsInROI -- The bounding box around the detected phantom is: X=({1} - {2}); Y=({3} - {4})'.format(os.path.basename(__file__), iLeft, iRight, iTop, iBottom))

        print('[{0}] :: GetPixelsInROI - getPixelsInROI -- Getting all pixel coordinates and corresponding PET values inside the ROI based on the configured radius and detected central pixel.'.format(os.path.basename(__file__)))
        for x in range(iLeft, iRight + 1):
            for y in range(iTop, iBottom + 1):
                dist = math.pow(self.Phantom_center_Column - x, 2.0) + math.pow(self.Phantom_center_Row - y, 2.0)
                if dist <= math.pow(radRoiPixels, 2.0):
                    self.pixelDictionary[x, y] = (self.dcmData.pixel_array[y, x] * self.dcmData.RescaleSlope)


    def getPixelValues(self):
        print('[{0}] :: GetPixelsInROI - getPixelValues -- Getting corresponding PET values of included ROI pixels.'.format(os.path.basename(__file__)))
        for key in self.pixelDictionary:
            self.pixelDictionary[key] = (self.dcmData.pixel_array[key[1], key[0]] * self.dcmData.RescaleSlope)


    def getRoiQuadrants(self):
        print('[{0}] :: GetPixelsInROI - getRoiQuadrants -- Sorting the included ROI pixels per quadrant.'.format(os.path.basename(__file__)))
        # Define central pixel
        x0 = self.Phantom_center_Column
        y0 = self.Phantom_center_Row

        # Create self.pixelDictionaryionary for each quadrant
        self.pixelDictionaryQ1 = {}
        self.pixelDictionaryQ2 = {}
        self.pixelDictionaryQ3 = {}
        self.pixelDictionaryQ4 = {}

        for key in self.pixelDictionary:
            if key[0] >= x0:
                if key[1] < y0:
                    self.pixelDictionaryQ1[key] = self.pixelDictionary[key]
                if key[1] >= y0:
                    self.pixelDictionaryQ2[key] = self.pixelDictionary[key]
            if key[0] <= x0:
                if key[1] < y0:
                    self.pixelDictionaryQ4[key] = self.pixelDictionary[key]
                if key[1] >= y0:
                    self.pixelDictionaryQ3[key] = self.pixelDictionary[key]


    def getROIMean(self):
        # get mean value of ROI
        self.ROImean = np.mean([self.pixelDictionary[key] for key in self.pixelDictionary])

        # Scale from Bq/ml to kBq/ml
        if self.dcmData.Units == 'BQML':
            self.ROImean = self.ROImean / 1000
            print('[{0}] :: GetPixelsInROI - getROIMean -- The mean value of the ROI (kBq/ml): {1}'.format(os.path.basename(__file__), self.ROImean))
        else:
            print('[{0}] :: GetPixelsInROI - getROIMean -- ERROR: The dicom data is stored in unexpected units!'.format(os.path.basename(__file__)))


    def getROIsd(self):
        # get standard deviation of ROI
        self.ROIsd = np.std([self.pixelDictionary[key] for key in self.pixelDictionary])

        # Scale from Bq/ml to kBq/ml
        if self.dcmData.Units == 'BQML':
            self.ROIsd = self.ROIsd / 1000
            print('[{0}] :: GetPixelsInROI - getROIsd -- The standard deviation of the ROI (kBq/ml): {1}'.format(os.path.basename(__file__), self.ROIsd))
        else:
            print('[{0}] :: GetPixelsInROI - getROIsd -- ERROR: The dicom data is stored in unexpected units!'.format(os.path.basename(__file__)))


    def getQROIMean(self):
        # get mean values of ROI
        self.ROImeanQ1 = np.mean([self.pixelDictionaryQ1[key] for key in self.pixelDictionaryQ1])
        self.ROImeanQ2 = np.mean([self.pixelDictionaryQ2[key] for key in self.pixelDictionaryQ2])
        self.ROImeanQ3 = np.mean([self.pixelDictionaryQ3[key] for key in self.pixelDictionaryQ3])
        self.ROImeanQ4 = np.mean([self.pixelDictionaryQ4[key] for key in self.pixelDictionaryQ4])

        # Scale from Bq/ml to kBq/ml
        if self.dcmData.Units == 'BQML':
            self.ROImeanQ1 = self.ROImeanQ1 / 1000
            self.ROImeanQ2 = self.ROImeanQ2 / 1000
            self.ROImeanQ3 = self.ROImeanQ3 / 1000
            self.ROImeanQ4 = self.ROImeanQ4 / 1000
            print('[{0}] :: GetPixelsInROI - getQROIMean -- The mean values of the ROI quadrants (kBq/ml): Q1={1}; Q2={2}; Q3={3}; Q4={4}'.format(os.path.basename(__file__), self.ROImeanQ1, self.ROImeanQ2, self.ROImeanQ3, self.ROImeanQ4))
        else:
            print('[{0}] :: GetPixelsInROI - getQROIMean -- ERROR: The dicom data is stored in unexpected units!'.format(os.path.basename(__file__)))


    def getQROIsd(self):
        # get standard deviation of ROI
        self.ROIsdQ1 = np.std([self.pixelDictionaryQ1[key] for key in self.pixelDictionaryQ1])
        self.ROIsdQ2 = np.std([self.pixelDictionaryQ2[key] for key in self.pixelDictionaryQ2])
        self.ROIsdQ3 = np.std([self.pixelDictionaryQ3[key] for key in self.pixelDictionaryQ3])
        self.ROIsdQ4 = np.std([self.pixelDictionaryQ4[key] for key in self.pixelDictionaryQ4])

        # Scale from Bq/ml to kBq/ml
        if self.dcmData.Units == 'BQML':
            self.ROIsdQ1 = self.ROIsdQ1 / 1000
            self.ROIsdQ2 = self.ROIsdQ2 / 1000
            self.ROIsdQ3 = self.ROIsdQ3 / 1000
            self.ROIsdQ4 = self.ROIsdQ4 / 1000
            print('[{0}] :: GetPixelsInROI - getQROIsd -- The standard deviation of the ROI quadrants (kBq/ml): Q1={1}; Q2={2}; Q3={3}; Q4={4}'.format(os.path.basename(__file__), self.ROIsdQ1, self.ROIsdQ2, self.ROIsdQ3, self.ROIsdQ4))
        else:
            print('[{0}] :: GetPixelsInROI - getQROIsd -- ERROR: The dicom data is stored in unexpected units!'.format(os.path.basename(__file__)))



class PerformAnalysis:

    def __init__(self, resultDictionary):
        self.results = resultDictionary


    def calculateUniformity(self, quadrant):
        listMean = []
        listSd = []

        for slc in self.results:
            listMean.append(self.results[slc][quadrant]['mean'])
            listSd.append(self.results[slc][quadrant]['std'])

        if min(listMean) > 0:
            self.Quniformity = np.mean(listSd)/np.mean(listMean)
            self.QuniformityPercentage = round((self.Quniformity * 100), 4)
            print('[{0}] :: PerformAnalysis - calculateUniformity -- The calculated uniformity for {1} is: {2}%'.format(os.path.basename(__file__), quadrant, self.QuniformityPercentage))

        else:
            print('[{0}] :: PerformAnalysis - calculateUniformity -- ERROR: Cannot calculate uniformity for {1} because not all mean ROI values are above zero!'.format(os.path.basename(__file__), quadrant))


    def calculateInSliceUniformity(self, slcKey):
        print('[{0}] :: PerformAnalysis - calculateInSliceUniformity -- Calculating in-slice uniformity for slice {1}.'.format(os.path.basename(__file__), slice))

        listMean = []
        listSd = []
        
        for quadrant in self.results[slcKey]:
            # skip global mean and SD
            if 'global' not in quadrant:
                listMean.append(self.results[slcKey][quadrant]['mean'])
                listSd.append(self.results[slcKey][quadrant]['std'])

        if min(listMean) > 0:
            sliceUniformity = np.mean(listSd)/np.mean(listMean)
            self.inSliceUniformity = round((sliceUniformity * 100), 4)
            print('[{0}] :: PerformAnalysis - calculateInSliceUniformity -- The calculated in-slice uniformity of slice {1} is: {2}%'.format(os.path.basename(__file__), slice, self.inSliceUniformity))

        else:
            print('[{0}] :: PerformAnalysis - calculateInSliceUniformity -- ERROR: Cannot calculate in-slice uniformity because not all mean ROI values are above zero!'.format(os.path.basename(__file__)))


    def calculateCalibrationDeviation(self, quadrant, corrConcentration):
        listMean = []

        for slc in self.results:
            listMean.append(self.results[slc][quadrant]['mean'])
        if min(listMean) > 0:
            self.calibrationDeviation = (np.mean(listMean) - corrConcentration) / corrConcentration
            self.calibrationDeviationPercentage = round((self.calibrationDeviation * 100), 4)
            print('[{0}] :: PerformAnalysis - calculateCalibrationDeviation -- The calculated deviation from calibration is: {1}%'.format(os.path.basename(__file__), self.calibrationDeviationPercentage))
        else:
            print('[{0}] :: PerformAnalysis - calculateCalibrationDeviation -- ERROR: Cannot calculate the deviation from calibration because not all mean ROI values are above zero.'.format(os.path.basename(__file__)))


class PhantomCalibration:

    def __init__(self, params):
        self.dateOfCalibration = params["calibration_date"]
        self.initialActivity = float(params["ref_activity"])
        self.calibrationCorrectionFactor = float(params["corrFactor_calibration"])
        self.halfLife = float(params["half_life"])
        self.volumePhantom = float(params["phantom_volume"])
        self.phantomName = params["phantom_name"]


    def getAcquisitionDateTime(self, file):
        dcmFile = dicom.read_file(file)

        self.acquisitionDate = dcmFile.AcquisitionDate
        self.acquisitionTime = dcmFile.AcquisitionTime[0:6]
        #self.acquisitionTime = dcmFile.AcquisitionDateTime[8:]

        print('[{0}] :: PhantomCalibration - getAcquisitionDateTime -- The acquisition is performed on {1} at {2}'.format(os.path.basename(__file__), self.acquisitionDate, self.acquisitionTime))


    def calculateCorrectedActivity(self):
        calibrationDate = datetime.date(int(self.dateOfCalibration[:4]), int(self.dateOfCalibration[4:6]), int(self.dateOfCalibration[-2:]))
        acqDate = datetime.date(int(self.acquisitionDate[:4]), int(self.acquisitionDate[4:6]), int(self.acquisitionDate[-2:]))

        diff = acqDate - calibrationDate
        self.correctedActivity = round(self.calibrationCorrectionFactor * self.initialActivity * math.exp(-math.log(2) * diff.days / self.halfLife), 2)

        print('[{0}] :: PhantomCalibration - calculateCorrectedActivity -- The calculated corrected activity is: {1} MBq'.format(os.path.basename(__file__), self.correctedActivity))


    def calculateCorrectedConcentration(self):
        self.correctedConcentration = round(self.correctedActivity / self.volumePhantom * 1000, 2)

        print('[{0}] :: PhantomCalibration - calculateCorrectedConcentration -- The calculated corrected concentration is: {1} kBq/ml'.format(os.path.basename(__file__), self.correctedConcentration))



class VerificationImages:

    def __init__(self, dict, file, threshold, centralPixelX, centralPixelY):
        self.dictionary = dict
        self.dcmFile = file
        self.usedThreshold = threshold
        self.centralPixelColumn = centralPixelX
        self.centralPixelRow = centralPixelY


    def getRoiBoundaryPixels(self):
        # Determine the boundary pixels of the defined ROI
        roiBoundaryPixelsX = []
        roiBoundaryPixelsY = []
        for key in self.dictionary:
            tmpKey1 = key[0] + 1, key[1]
            tmpKey2 = key[0] - 1, key[1]
            tmpKey3 = key[0], key[1] + 1
            tmpKey4 = key[0], key[1] - 1

            if tmpKey1 in self.dictionary and tmpKey2 in self.dictionary and tmpKey3 in self.dictionary and tmpKey4 in self.dictionary:
                pass
            else:
                roiBoundaryPixelsX.append(key[0])
                roiBoundaryPixelsY.append(key[1])

        return roiBoundaryPixelsX, roiBoundaryPixelsY


    def createQuadrantpng(self, plotName, ROIvalues, ROI, slcKey, dcmInstancNumber):
        print('[{0}] :: VerificationImages - createQuadrantpng -- Creating verification image of ROI quadrants'.format(os.path.basename(__file__)))

        # dictionary with axis for each quadrant
        axq = { 'Q1': None, 'Q2': None, 'Q3': None, 'Q4': None }
        # define plot
        #f, ((ax4, ax1), (ax3, ax2)) = plt.subplots(2, 2, figsize=(7, 7))
        f, ((axq['Q4'], axq['Q1']), (axq['Q3'], axq['Q2'])) = plt.subplots(2, 2, figsize=(7, 7))
        titleSub = 'This slice (' + str(slcKey) + ') corresponds to DICOM Instance Number ' + str(dcmInstancNumber) + '.'
        f.suptitle(titleSub, color='red')

        # read dicom file
        ds = dicom.dcmread(self.dcmFile)

        # Get reference image
        refImage_2d = ds.pixel_array.astype(int)
        
        # Need to access pixel coordinates
        # These available as dictionary with coordinates (tuples) as keys
        # Would be nice to have these as a dictionary in the first place
        # but this is all over the code, so leave it for now.
        pixelDicts = { 'Q1': ROI.pixelDictionaryQ1, 'Q2': ROI.pixelDictionaryQ2, 'Q3': ROI.pixelDictionaryQ3, 'Q4': ROI.pixelDictionaryQ4 }

        # Loop over quadrants
        for q in axq:
            # SUB FIGURE
            xList = []
            yList = []
            # here key is actually a tuple with coordinates of pixel
            for key in pixelDicts[q]:
                xList.append(key[0])
                yList.append(key[1])
            axq[q].imshow(refImage_2d, origin='upper', cmap='gray', clim=None)
            axq[q].scatter(xList, yList, color='red', s=1, alpha=0.3)
            axq[q].text(10, 15, q, color='red', size=12)
            tmpText1 = 'Mean [kBq/ml]: ' + str(round(ROIvalues[slcKey][q]['mean'], 4))
            tmpText2 = 'SD   [kBq/ml]: ' + str(round(ROIvalues[slcKey][q]['std'] , 4))
            axq[q].text(10, 125, tmpText1, color='white', size=9)
            axq[q].text(10, 135, tmpText2, color='white', size=9)
            axq[q].axis('off')
        
        plt.savefig(plotName)
        #plt.close()


    def createROIpng(self, plotName, figName1, figName2, figName3, Key, Value, ROImean, ROIsd):
        print('[{0}] :: VerificationImages - createROIpng -- Creating verification image of ROI'.format(os.path.basename(__file__)))

        # define plot
        f, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(8, 4))
        titleSub = 'This slice (' + str(Key) + ') corresponds to DICOM Instance Number ' + str(Value) + '.'
        f.suptitle(titleSub, color='red')

        # read dicom file
        ds = dicom.dcmread(self.dcmFile)

        # Get reference image
        refImage_2d = ds.pixel_array.astype(int)

        # SUB FIGURE 1
        ax1.imshow(refImage_2d, origin='upper', cmap='gray', clim=None)
        ax1.axis('off')
        ax1.set_title(figName1, size=11, y=1.08)

        # SUB FIGURE 2
        xList = []
        yList = []
        for key in self.dictionary:
            xList.append(key[0])
            yList.append(key[1])
        ax2.imshow(refImage_2d, origin='upper', cmap='gray', clim=None)
        ax2.scatter(xList, yList, color='red', s=1, alpha=0.3)
        ax2.scatter(self.centralPixelColumn, self.centralPixelRow, color='blue', s=1)
        ax2.text(10, 10, 'Detected central pixel', color='blue', size=8)
        ax2.text(10, 25, 'Detected ROI', color='red', size=8)
        ax2.text(10, 130, 'ROI mean: {:.2f} kBq/ml'.format(ROImean), color='white', size=8)
        ax2.text(10, 120, 'ROI sd: {:.3f} kBq/ml'.format(ROIsd), color='white', size=8)
        ax2.axis('off')
        ax2.set_title(figName2, size=11, y=1.08)

        # SUB FIGURE 3
        roiX, roiY = self.getRoiBoundaryPixels()
        ax3.imshow(refImage_2d, origin='upper', cmap='gray', clim=None)
        ax3.scatter(self.centralPixelColumn, self.centralPixelRow, color='blue', s=1)
        ax3.scatter(roiX, roiY, color='red', s=1, alpha=0.3)
        ax3.text(10, 10, 'Detected central pixel', color='blue', size=8)
        ax3.text(10, 25, 'Detected ROI (contour)', color='red', size=8)
        ax3.axis('off')
        ax3.set_title(figName3, size=11, y=1.08)

        plt.savefig(plotName)
        #plt.close()


    def createMIP(self, data, slices, plotName, figTitle):
        # define number of additional rows in image
        additionalRows = 5

        # define size of mip
        dcmFile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)
        sizeX = dcmFile.Columns
        sizeY = dcmFile.NumberOfSlices + (additionalRows*3)

        # define empty MIP image
        MIP = [[0 for x in range(sizeX)] for y in range(sizeY)]

        # loop over dicom files
        dcmFilesList = data.series_filelist[0]
        for file in dcmFilesList:

            # get pixel data from instance
            instance = dicom.read_file(file)
            pixels = instance.pixel_array.astype(int)
            mipRow = np.amax(pixels, axis=0)

            # copy mipRow values to MIP image
            for i in range(len(mipRow)):
                MIP[(instance.InstanceNumber + additionalRows)][i] = mipRow[i]

        # define plot
        fig, ax = plt.subplots()
        ax.imshow(MIP, origin='lower', cmap='gray', clim=None)

        # add slice numbers used for analysis
        for key in slices:
            plt.hlines(slices[key]+additionalRows, 30, sizeX-30, colors='red')
            plt.text(22, slices[key]+additionalRows-1, str(key), color='red')
        plt.text(10, 55, 'analyzed slices:', color='red')

        # add axis, labels, title, etc..
        ax.set(xlabel='X [pixels]', ylabel='DICOM Instance Number')
        plt.xticks([])
        plt.yticks((5, 20, 35, 50), ('0', '15', '30', '45'))
        ax.set_title(figTitle, size=10, y=1.08)

        # save image
        plt.savefig(plotName)
        #plt.close()
