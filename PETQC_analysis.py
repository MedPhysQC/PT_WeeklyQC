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
# PyWAD is open-source software and consists of a set of modules written in python for the WAD-Software medical physics quality control software.
# The WAD Software can be found on https://github.com/wadqc
#
# The pywad package includes modules for the automated analysis of QC images for various imaging modalities.
# PyWAD has been originaly initiated by Dennis Dickerscheid (AZN), Arnold Schilham (UMCU), Rob van Rooij (UMCU) and Tim de Wit (AMC)


"""

Dependencies:

    PYTHON
    - numpy
    - pydicom / dicom
    - datetime
    - ast


Version control:
    20181128: first version
    20190227: [JK] add deviation from calibration per quadrant to results
    20190228: [JK] reverse output images (was upside down), change slc/quadant/metric indexing
    20190306: [JK] fix issue with text on quadrant images
    20220916: [MY] add auto slice selection
    20220928: [MY] bug fixes auto slice selection

The initial version was created by Stijn van de Schoot at the Amsterdam UMC, location VUmc

"""

__version__ = '20220928'
__author__ = 'Stijn van de Schoot, Maqsood Yaqub'
__location__ = 'Amsterdam UMC'


#import sys
#import math
import os
#import numpy as np
#import datetime
import ast
import math # 20220916 MY added for auto slice routine
import numpy # 20220916 MY added for auto slice routine

# Import WAD QC modules
from wad_qc.module import pyWADinput
from wad_qc.modulelibs import wadwrapper_lib


# Import required module(s)
try:
    import PETQC_lib
except ImportError:
    print('[{0}] :: Could not import PETQC_lib module.'.format(os.path.basename(__file__)))

# Import DICOM
try:
    import pydicom as dicom
except ImportError:
    import dicom




def logTag():
    return "[{0}] ::".format(os.path.basename(__file__))

def logTagInitialisation():
    return "[****]"

def startLogging():
    print('===========================================================================================')
    print('{0} The WAD QC analysis module is written by {1} ({2})'.format(logTagInitialisation(), __author__, __location__))
    print('{0} The current version used is: {1}'.format(logTagInitialisation(), __version__))
    print('\n')
    print('{0} This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;'.format(logTagInitialisation()))
    print('{0} without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'.format(logTagInitialisation()))
    print('{0} See the GNU General Public License for more details.'.format(logTagInitialisation()))
    print('\n')
    print('{0} Included module(s):'.format(logTagInitialisation()))
    print('{0}      - {1} (version: {2}; author: {3})'.format(logTagInitialisation(), PETQC_lib.__name__, PETQC_lib.__version__, PETQC_lib.__author__))
    print('===========================================================================================')
    print('\n')




def getAcquistionDateTime(data, results):

    """
    Read acquisition data and time from dicomheaders and write to IQC database

    Workflow:
        1. Read only header and get date and time
    """

    ## 1. read only header and get date and time

    # Add info to logfile
    print('{} getAcquisitionDateTime -- Reading DICOM header and get date and time of acquisition.'.format(logTag()))

    # Get DICOM info
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)
    dt = wadwrapper_lib.acqdatetime_series(dcmInfile)

    # Add data to results
    print('{} getAcquisitionDateTime -- Adding date and time of acquisition to results.'.format(logTag()))
    print('\n')
    results.addDateTime('Acquisition date and time', dt)



def calculatePhantomActivity(data, results, action):

    """
    Calculate current activity of phantom based on calibration data from the module configuration file and the scan date and time from the dicom file.

    Workflow:
        1. Initiate class and get parameter info
        2. Get acquisition date and time
        3. Calculate corrected activity
        4. Calculate corrected concentration
    """

    try:
        params = action['params']
    except KeyError:
        print('{} calculatePhantomActivity -- Could not read parameters from config file.'.format(logTag()))
        params = {}

    # Verify content of parameter dictionary
    if len(params) == 0:
        print('{} calculatePhantomActivity -- ERROR: no parameters found!'.format(logTag()))

    # Define output dictionary
    phantomData = {}

    ## 1, Initiate class
    print('{} calculatePhantomActivity -- Initiating class.'.format(logTag()))
    Correction = PETQC_lib.PhantomCalibration(params)
    results.addString("Phantom ID", Correction.phantomName)
    results.addString("Phantom calibration date", Correction.dateOfCalibration)
    results.addFloat("Phantom reference activity", Correction.initialActivity)
    phantomData["Phantom ID"] = Correction.phantomName

    ## 2. Get acquisition date and time
    print('{} calculatePhantomActivity -- Getting acquisition date and time.'.format(logTag()))
    Correction.getAcquisitionDateTime(data.series_filelist[0][0])

    ## 3. Calculate corrected activity
    print('{} calculatePhantomActivity -- Calculating corrected activity.'.format(logTag()))
    Correction.calculateCorrectedActivity()
    results.addFloat("Corrected Activity (MBq)", Correction.correctedActivity)
    phantomData["Corrected Activity (MBq)"] = Correction.correctedActivity

    ## 4. Calculate corrected concentration
    print('{} calculatePhantomActivity -- Calculating corrected concentration.'.format(logTag()))
    Correction.calculateCorrectedConcentration()
    results.addFloat("Corrected Concentration (kBq/ml)", Correction.correctedConcentration)
    phantomData["Corrected Concentration (kBq/ml)"] = Correction.correctedConcentration

    print('\n')
    return phantomData

def getSlicesFromLocationPhantomInImage(listoffiles):
    # 20220916 MY auto slice locating routine

    # checking filenames
    files = []
    for fname in listoffiles:
        files.append(dicom.dcmread(fname)) # future need to catch non dicom file error
    slices = []
    skipcount = 0
    for f in files:
        if hasattr(f, 'SliceLocation'):
            slices.append(f)
        else:
            skipcount = skipcount + 1

    slices = sorted(slices, key=lambda s: s.SliceLocation)

    # convert to filenames
    filessorted = []
    for file in slices:
        filessorted.append(file.filename)

    minvalue = float ('1000000')
    maxvalue = float ('-1000000')

    # locate slices with phantom
    for file in filessorted:
        dcmFile = dicom.read_file(file)
        data = dcmFile.pixel_array
        meanvalue= data.mean()*float(dcmFile.RescaleSlope)
        if(math.isfinite(meanvalue)):
            if  meanvalue > maxvalue:
                maxvalue = meanvalue
            if  meanvalue < minvalue:
                minvalue = meanvalue

    threshold = maxvalue * 0.2;
    idealInstanceNumbers = []; # number of slices with phantom
    for file in filessorted:
        dcmFile = dicom.read_file(file)
        data = dcmFile.pixel_array
        meanvalue=data.mean()*float(dcmFile.RescaleSlope)
        if(math.isfinite(meanvalue)):
            if  meanvalue > threshold:
                idealInstanceNumbers.append(int(dcmFile.InstanceNumber))

    #print(idealInstanceNumbers)
    minslice = numpy.min(idealInstanceNumbers)
    maxslice = numpy.max(idealInstanceNumbers)

    stepsize = int(math.floor((maxslice - minslice)/4));
    sliceA = minslice+stepsize;
    sliceB = sliceA+stepsize;
    sliceC = sliceB+stepsize;

    return [sliceA, sliceB, sliceC]

def petQcAnalysis(data, results, action, phantomData):

    """
    Perform PET QC analysis

    Workflow:
        1. Get included DICOM Instance Numbers
        2. For each included slice:
            - get ROI and ROI quadrant statistics
            - create image of ROI and ROI quadrants
        3. Create maximum intensity projection image
        4. Calculate uniformity (entire ROI and per ROI quadrant) and deviation from calibration
    """

    try:
        params = action['params']
    except KeyError:
        print('{} petQcAnalysis -- ERROR: Could not read parameters from config file.'.format(logTag()))
        params = {}

    # Verify content of parameter dictionary
    if len(params) == 0:
        print('{} petQcAnalysis -- ERROR: no parameters found!'.format(logTag()))


    ## 1. Get list of included DICOM Instance Numbers
    print('{} petQcAnalysis -- Getting included DICOM Instance Numbers from config.'.format(logTag()))
    # 20220928 switch for auto slice detection
    if params["included_slices"] == 'autodetect':
        # autodetect slices
        configuredDcmIN = getSlicesFromLocationPhantomInImage(data.series_filelist[0])
        stringDcmIN = '{} petQcAnalysis -- Auto-detected slices: ' + str(configuredDcmIN)
    else:
        # use configured slices
        configuredDcmIN = ast.literal_eval(params["included_slices"])
        stringDcmIN = '{} petQcAnalysis -- Configured slices: ' + str(configuredDcmIN)
    print(stringDcmIN.format(logTag()))

    print('{} petQcAnalysis -- Assigning unique letters to included DICOM Instance Numbers from config.'.format(logTag()))
    dcmInstanceNumbers = {}
    index = PETQC_lib.IndexConverter()
    for i, dcmIN in enumerate(configuredDcmIN):
        tmpKey = index.indexToLetter(i)
        dcmInstanceNumbers[tmpKey] = dcmIN

        # Add conversion to results
        strOutput = 'Slice ' + str(tmpKey)
        results.addFloat(strOutput, dcmIN)


    ## 2. For each included DICOM Instance Number, get ROI statistics
    ROIvalues = {}

    # loop over requested (configured) slices
    for slcKey in dcmInstanceNumbers:
        print('{0} petQcAnalysis -- Starting analysis for slice {1} corresponding to DICOM Instance Number {2}.'.format(logTag(), slcKey, dcmInstanceNumbers[slcKey]))

        # loop over dicom files
        dcmFilesList = data.series_filelist[0]
        for file in dcmFilesList:
            dcmFile = dicom.read_file(file, stop_before_pixels=True)

            if dcmInstanceNumbers[slcKey] == dcmFile.InstanceNumber:
                print('{0} petQcAnalysis -- Dicom file for slice {1} corresponding to DICOM Instance Number {2} is found.'.format(logTag(), slcKey, dcmInstanceNumbers[slcKey]))

                # Initiate class
                ROI = PETQC_lib.ROI(params, file)

                # Determine the central pixel of the phantom
                ROI.getCentralPixel()

                # Determine all pixels in the ROI based on the central pixel and the ROI radius and collect corresponding PET values
                ROI.getPixelsInROI()

                # Determine mean and sd value of ROI
                ROI.getROIMean()
                outputString = 'Slice ' + str(slcKey) + ' - ROI mean'
                results.addFloat(outputString, ROI.ROImean)

                ROI.getROIsd()
                outputString = 'Slice ' + str(slcKey) + ' - ROI sd'
                results.addFloat(outputString, ROI.ROIsd)

                # Store values in dictionary
                ROIvalues[slcKey] = {}
                ROIvalues[slcKey]['global'] = {}
                ROIvalues[slcKey]['global']['mean'] = ROI.ROImean
                ROIvalues[slcKey]['global']['std']  = ROI.ROIsd

                # Create image of analysed ROI
                images = PETQC_lib.VerificationImages(ROI.pixelDictionary, file, ROI.threshold,
                                                      ROI.Phantom_center_Column, ROI.Phantom_center_Row)
                plotName = 'PETQC_slice_' + str(slcKey) + '.png'
                titleFig1 = 'Original slice'
                titleFig2 = 'ROI'
                titleFig3 = 'ROI contour'
                images.createROIpng(plotName, titleFig1, titleFig2, titleFig3, slcKey, dcmInstanceNumbers[slcKey], ROI.ROImean, ROI.ROIsd)
                resultName = 'Slice ' + str(slcKey) + ' - ROI visualisation'
                results.addObject(resultName, plotName)


                ## ROI QUADRANTS
                # Split pixel dictionary per quadrant
                # TODO: ROI quadrants could be ROI objects themselves,
                #       and then we could have a dict with ROIs for indexed access
                ROI.getRoiQuadrants()

                # Determine mean and sd value of ROI quadrants
                ROI.getQROIMean()
                ROI.getQROIsd()

                # Store values in dictionary
                ROIvalues[slcKey]['Q1'] = {}
                ROIvalues[slcKey]['Q1']['mean'] = ROI.ROImeanQ1
                ROIvalues[slcKey]['Q1']['std'] = ROI.ROIsdQ1
                ROIvalues[slcKey]['Q2'] = {}
                ROIvalues[slcKey]['Q2']['mean'] = ROI.ROImeanQ2
                ROIvalues[slcKey]['Q2']['std'] = ROI.ROIsdQ2
                ROIvalues[slcKey]['Q3'] = {}
                ROIvalues[slcKey]['Q3']['mean'] = ROI.ROImeanQ3
                ROIvalues[slcKey]['Q3']['std'] = ROI.ROIsdQ3
                ROIvalues[slcKey]['Q4'] = {}
                ROIvalues[slcKey]['Q4']['mean'] = ROI.ROImeanQ4
                ROIvalues[slcKey]['Q4']['std'] = ROI.ROIsdQ4

                # Calculate in-slice uniformity
                Statistics = PETQC_lib.PerformAnalysis(ROIvalues)
                Statistics.calculateInSliceUniformity(slcKey)
                resultsName = 'Slice ' + str(slcKey) + ' - in-slice uniformity'
                results.addFloat(resultsName, Statistics.inSliceUniformity)

                # Create image of ROI quadrants
                plotName = 'PETQC_slice_' + str(slcKey) + '_Quadrants.png'
                images.createQuadrantpng(plotName, ROIvalues, ROI, slcKey, dcmInstanceNumbers[slcKey])
                resultName = 'Slice ' + str(slcKey) + ' - ROI Quadrants visualisation'
                results.addObject(resultName, plotName)


    # add input data to results
    results.addFloat('ROI radius', ROI.radRoiCM)
    results.addFloat('threshold value', ROI.threshold)

    print('{0} :: petQcAnalysis -- The extracted ROI data are as follows: '.format(logTag()))
    for slc in ROIvalues:
        for quadrant in ROIvalues[slc]:
            for metric in ROIvalues[slc][quadrant]:
                print('{0} :: petQcAnalysis -- slice {1} (#{2}) {3} {4} = {5}'.format(logTag(), slc, dcmInstanceNumbers[slc], quadrant, metric, ROIvalues[slc][quadrant][metric]))


    ## 3. Create MIP image
    print('{} petQcAnalysis -- Creating Maximum Intensity Projection (MIP) image.'.format(logTag()))
    plotMIPname = 'PETQC_MIP.png'
    mipTitle = 'Maximum intensity projection (MIP)'
    images.createMIP(data, dcmInstanceNumbers, plotMIPname, mipTitle)
    results.addObject("MIP image", plotMIPname)


    ## 4. Calculate uniformity and deviation from calibration
    print('{} petQcAnalysis -- Calculating uniformity.'.format(logTag()))
    Statistics = PETQC_lib.PerformAnalysis(ROIvalues)
    quadrantList = ['global', 'Q1', 'Q2', 'Q3', 'Q4']
    for q in quadrantList:
        print('{0} petQcAnalysis -- Start uniformity calculation for ROI: {1}.'.format(logTag(), q))
        Statistics.calculateUniformity(q)
        results.addFloat(q + ' uniformity', Statistics.QuniformityPercentage)

        print('{0} petQcAnalysis -- Calculating deviation from calibration.'.format(logTag(), q))
        Statistics.calculateCalibrationDeviation(q, phantomData["Corrected Concentration (kBq/ml)"])
        results.addFloat(q + ' deviation from calibration', Statistics.calibrationDeviationPercentage)



def dataVerification(data, actions):

    # read parameters
    try:
        action = actions['dataVerification']
        params = action['params']
    except KeyError:
        print('{} dataVerification -- ERROR: Could not read parameters from config file.'.format(logTag()))
        params = {}

    # define initial result
    result = True

    # read the header of the first DICOM file
    dcmFile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    # Check if pixels are equally spaced in both dimensions
    if dcmFile.PixelSpacing[0] != dcmFile.PixelSpacing[1]:
        print('{} dataVerification -- Incorrect pixel spacing.'.format(logTag()))
        result = False

    # Check if correct PET units are used
    if dcmFile.Units != params['units']:
        print('{} dataVerification -- Incorrect PET units.'.format(logTag()))
        result = False
    results.addString("Pixel units (DICOM)", dcmFile.Units)

    # Check image type
    if dcmFile.ImageType[0] != params['image_type'][0] or dcmFile.ImageType[1] != params['image_type'][1]:
        print('{} dataVerification -- Unknown Image Type.'.format(logTag()))
        result = False
    strResult = params['image_type'][0] + '/' + params['image_type'][1]
    results.addString("Image type (DICOM)", strResult)

    # Check slice thickness
    if dcmFile.SliceThickness != float(params['slice_thickness']):
        print('{} dataVerification -- Unknown Slice Thickness.'.format(logTag()))
        result = False
    results.addFloat("Slice thickness (DICOM)", dcmFile.SliceThickness)

    # Check protocol name
    if dcmFile.ProtocolName != params['protocol_name']:
        print('{} dataVerification -- Incorrect Protocol Name.'.format(logTag()))
        result = False
    results.addString("Protocol name (DICOM)", dcmFile.ProtocolName)

    return result



if __name__=="__main__":

    # start logging
    startLogging()

    # initialisation of WAD module
    data, results, config = pyWADinput()

    # get acquisition time and date
    for name, action in config['actions'].items():
        if name == 'acqDateTime':
            getAcquistionDateTime(data, results)

    # verify incoming DICOM data
    if dataVerification(data, config['actions']):

        # calculate the phantom activity
        for name, action in config['actions'].items():
            if name == 'phantomActivity':
                phantomData = calculatePhantomActivity(data, results, action)

        # perform the actual PET QC analysis
        for name, action in config['actions'].items():
            if name == 'petQC_analysis':
                petQcAnalysis(data, results, action, phantomData)

    else:
        print('{} Main -- ERROR: data does not fulfill all requirements!'.format(logTag()))

    # write results to database
    results.write()
