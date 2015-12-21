from __future__ import division

import re
import logging
from glob import glob

import numpy as np
import h5py

#from pbpy.align.ZScore import ZScoreCalculator

log = logging.getLogger(__name__)

__version__ = '1.1'

COMPUTED_FIELDS = ['StartTime', 'PulseWidth', 'pkmid', 'pkmax', 'IPD',
                   'Light']

FIELD_DEPENDENCIES = {'StartTime': (['StartFrame'], []),
                      'PulseWidth': (['WidthInFrames'], ['WidthInFrames']),
                      'pkmid': (['MidSignal', 'Channel'], []),
                      'pkmax': (['MaxSignal', 'Channel'], []),
                      'pkmean': (['MeanSignal', 'Channel'], []),
                      'IPD': (['StartFrame', 'WidthInFrames'], ['PreBaseFrames']),
                      'Light': (['WidthInFrames', 'MeanSignal', 'Channel'], [])
                      }

PLS_REGION_TYPES = ['HQRegion']
PLS_REGION_DESC = ['Denotes highest quality wells']
PLS_REGION_SRC = ['SMRTpipe Filter Module']


class MissingMetricError(Exception):
    pass


class PlsRegion(object):

    """Defines an entry in the Regions table of a pls.h5 file."""

    TABLE_COLUMNS = ["HoleNumber", "TypeIndex", "Start", "End", "Score"]

    def __init__(self, **kwargs):
        self._data = {}
        for key in PlsRegion.TABLE_COLUMNS:
            self[key] = kwargs.get(key, None)

    def __setitem__(self, key, value):
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __len__(self):
        return self["End"] - self["Start"]

    def toTableRow(self):
        """Converts this Region into a valid entry in the regions
        table of a pls.h5 file."""
        return np.array([self._data[c] for c in PlsRegion.TABLE_COLUMNS], dtype=np.int32)

    def __str__(self):
        return "\n".join(["%s: %s" % (k, str(v)) for k, v in self._data.iteritems()])


def writeRegionsTable(regions, fileName, types=PLS_REGION_TYPES,
                      descriptions=PLS_REGION_DESC, sources=PLS_REGION_SRC):
    """Writes out a pls.h5 file containing a regions table defined
    by the arguments to this function."""

    outFile = h5py.File(fileName, 'w')

    # Add Version to root level of rgn.h5 File
    outFile.attrs['Version'] = __version__
    log.debug("Writing VERISON {x}".format(x=__version__))

    shape = (max(1, len(regions)), len(PlsRegion.TABLE_COLUMNS))
    pd = outFile.create_group("PulseData")
    regionTable = pd.create_dataset("Regions", shape, np.int32, maxshape=(None, shape[1]))

    rows = np.zeros(shape=shape, dtype=np.int32)
    for i, row in enumerate([region.toTableRow() for region in regions]):
        rows[i] = row
    regionTable[:] = rows

    local_types = [f for f in types]
    local_descriptions = [f for f in descriptions]
    local_sources = [f for f in sources]
    local_types = [str(f) for f in local_types]
    local_descriptions = [str(f) for f in local_descriptions]
    local_sources = [str(f) for f in local_sources]
    regionTable.attrs["ColumnNames"] = np.array(
        PlsRegion.TABLE_COLUMNS, dtype=h5py.new_vlen(str))
    regionTable.attrs["RegionTypes"] = np.array(local_types,
                                                   dtype=h5py.new_vlen(str))
    regionTable.attrs["RegionDescriptions"] = np.array(
        local_descriptions, dtype=h5py.new_vlen(str))
    regionTable.attrs["RegionSources"] = np.array(
        local_sources, dtype=h5py.new_vlen(str))

    outFile.close()


def regionTypes(h5FN):
    """Returns the list of region types from a pls file."""
    h5File = h5py.File(h5FN, 'r')
    if "/PulseData/Regions" in h5File and "RegionTypes" in h5File["/PulseData/Regions"].attrs:
        rTypes = h5File["/PulseData/Regions"].attrs["RegionTypes"]
    else:
        rTypes = []
    h5File.close()
    return rTypes


def regionDescriptions(h5FN):
    """Returns the list of region descriptions from a pls file."""
    h5File = h5py.File(h5FN, 'r')
    if "/PulseData/Regions" in h5File and "RegionDescriptions" in h5File["/PulseData/Regions"].attrs:
        rTypes = h5File["/PulseData/Regions"].attrs["RegionDescriptions"]
    else:
        rTypes = []
    h5File.close()
    return rTypes


def regionSources(h5FN):
    """Returns the list of region sources from a pls file."""
    h5File = h5py.File(h5FN, 'r')
    if "/PulseData/Regions" in h5File and "RegionSources" in h5File["/PulseData/Regions"].attrs:
        rTypes = h5File["/PulseData/Regions"].attrs["RegionSources"]
    else:
        rTypes = []
    h5File.close()
    return rTypes


def regionIterator(h5FN, rTypesToKeep=()):
    """
    Iterates over PlsRegion objects representing all regions of the given
    type. Type can be specified as a string (RegionType), an
    int (RegionTypeId), or a tuple of strings or ints.
    """

    if isinstance(rTypesToKeep, str) or isinstance(rTypesToKeep, int):
        rTypesToKeep = (rTypesToKeep, )

    h5File = h5py.File(h5FN, 'r')

    if "/PulseData/Regions" not in h5File:
        raise StopIteration

    rgnTable = h5File["/PulseData/Regions"]
    rgnTypes = rgnTable.attrs["RegionTypes"]

    if len(rTypesToKeep) > 0 and isinstance(rTypesToKeep[0], str):
        str2int = dict([(rType, i) for i, rType in enumerate(rgnTypes)])
        rTypesToKeep = [str2int[t] for t in rTypesToKeep if t in str2int]

    myRgnTable = rgnTable.value
    myRgnTable.dtype = [(col, np.int32) for col in PlsRegion.
                        TABLE_COLUMNS]
    if len(rTypesToKeep) > 0:
        for rType in rTypesToKeep:
            returnRows = myRgnTable[myRgnTable["TypeIndex"] == rType]
            for row in returnRows:
                yield PlsRegion(**dict(zip(PlsRegion.TABLE_COLUMNS, row)))  # passing the dict as kwargs
    else:
        for row in myRgnTable:
            yield PlsRegion(**dict(zip(PlsRegion.TABLE_COLUMNS, row[0])))  # passing the dict as kwargs

    h5File.close()


def regionMerge(rPaths, targetPath):
    """Given a list of paths to regions pls.h5 files,
    combines them into a single regions file."""

    regions = []
    for path in rPaths:
        rH5 = h5py.File(path, 'r')
        for row in rH5["/PulseData/Regions"]:
            regions.append(PlsRegion(HoleNumber=row[0],
                                     TypeIndex=row[1],
                                     Start=row[2],
                                     End=row[3],
                                     Score=row[4]))

    writeRegionsTable(regions, targetPath)


class MovieMetaData(object):

    def __init__(self, h5FN):  # hdf5FN: the *.mmd.h5 file name
        self.h5handle = None
        try:
            self.h5handle = h5py.File(h5FN, 'r')
        except IOError:
            raise IOError("Can not open HDF5 file: %s!!" % h5FN)
        self._metaData = {}
        self._fillMetaData()

    # The following two method fill the attribute from the metafile into an
    # internal dictionary
    def _fillOneAttribute(self, path, attr, data, attrDict):
        path = path.strip().split("/")
        if path != ['']:  # not leaves
            attrDict.setdefault(path[0], {})
            self._fillOneAttribute(
                "/".join(path[1:]), attr, data, attrDict[path[0]])
        else:
            attrDict["_" + attr] = data

    def _addAttribute(self, groupName, groupObj):
        for attr in groupObj.attrs:
            p = groupObj.name[1:]
            data = groupObj.attrs[attr]
            if not np.isscalar(data):
                data = data[0]
            self._fillOneAttribute(p, attr, data, self._metaData)

    def _fillMetaData(self):
        metaData = {}
        # TODO check if ScanData is in the file, if not, freak out
        scanData = self.h5handle["/ScanData"]
        scanData.visititems(self._addAttribute)
        # allow ChipMask access to gracefully fail (handles early Springfield
        # files)
        try:
            self._chipMask = self.h5handle["/ScanData/ChipArray/ChipMask"]
        except Exception, e:
            log.error(str(e))
            pass
        # print >>sys.stderr, "WARNING %s: couldn't access
        # /ScanData/ChipArray/ChipMask" % self.h5handle.filename
        # allow HolesPerLine access to gracefully fail (handles early
        # Springfield files)
        try:
            # avoid frequently lookup through pytables lib
            self.__HolesPerLine = self._HolesPerLine.values()[0]
        except Exception, e:
        # print >>sys.stderr, "WARNING %s: couldn't access HolesPerLine" %
        # self.h5handle.filename
            # for early Springfield data we need to dummy this
            self.__HolesPerLine = -1

            self.__HoleNumber = self.h5handle[
                "/PulseData/BaseCalls/ZMW/HoleNumber"][:]
            self.__HoleXY = self.h5handle["/PulseData/BaseCalls/ZMW/HoleXY"][:]
            self.__HoleXY2HoleNumber = dict(zip([(xy[0], xy[
                1]) for xy in self.__HoleXY], self.__HoleNumber))
            self.__HoleNumber2HoleXY = dict(zip(self.
                                                __HoleNumber, [(xy[0], xy[1]) for xy in self.__HoleXY]))
        # print >>sys.stderr, "WARNING %s: Calculated XY and hole number using
        # dictionary for SpringField" % self.h5handle.filename

    # Internal utility function, it returns key(path),value(data) pair for an
    # attribute name
    def _findAttrib(self, attr, attrD, path, collection=[]):
        if type(attrD) != type({}):
            return []
        if attr in attrD:
            path = path + "/" + attr
            return collection + [(path, attrD[attr])]
        else:
            rtn = []
            for a in attrD:
                p = path + "/" + a
                rtn += self._findAttrib(attr, attrD[a], p, collection)
            return rtn

    # Python's magic
    def __getattr__(self, attr):
        r = self._findAttrib(attr, self._metaData, "")
        if r == []:
            return None
        else:
            return dict(r)

    def XY2HoleNumber(self, x, y):
        if self.__HolesPerLine > 0:
            return self.__HolesPerLine * (x - 1) + (y - 1)
        else:
            return self.__HoleXY2HoleNumber[(x, y)]

    def holeNumber2XY(self, n):
        if self.__HolesPerLine > 0:
            x = int(n / self.__HolesPerLine) + 1
            y = int(n % self.__HolesPerLine) + 1
        else:
            (x, y) = self.__HoleNumber2HoleXY[n]
        return (x, y)

    # The following is to conform with the PacBio.Analysis.Data.MovieMetaZMWer
    @property
    def BaseMap(self):
        return self._BaseMap.values()[0]

    @property
    def Wavelength(self):
        w = self._Wavelength.items()
        w.sort()
        return [x[1] for x in w]

    # The following is to conform with the PacBio.Analysis.Data.MovieMetaZMWer
    @property
    def movieName(self):
        return self._MovieName.values()[0]

    @property
    def RunCode(self):
        return self._RunCode.values()[0]

    @property
    def InstrumentName(self):
        return self._InstrumentName.values()[0]

    @property
    def FrameRate(self):
        return self._FrameRate.values()[0]

    @property
    def NumFrames(self):
        return self._NumFrames.values()[0]

    @property
    def LaserOnFrame(self):
        return self._LaserOnFrame.values()[0]

    @property
    def HotStartFrame(self):
        return self._HotStartFrame.values()[0]

    @property
    def LaserIntensity(self):
        return self._LaserIntensity.values()[0]

    @property
    def CameraGain(self):
        return self._CameraGain.values()[0]

    @property
    def AduGain(self):
        return self._AduGain.values()[0]

    def __del__(self):
        if self.h5handle:
            self.h5handle.close()


class HDFZwmPulses(object):
    pass


class MoviePulseData(object):

    def __init__(self, h5FN, metadata=None, cachePulseData=True):  # hdf5FN: the *.pls.h5 file name, metadata: MovieMetaData instance
        self.h5FN = h5FN
        self.h5handle = None
        try:
            self.h5handle = h5py.File(h5FN, 'r')
        except IOError:
            raise IOError("Can not open HDF5 file: %s!!" % h5FN)
        if metadata is None:
            self._metaData = MovieMetaData(h5FN)
        else:
            self._metaData = metadata

        self.channel2base = dict(zip(range(4), str(self._metaData.BaseMap)))
        if 'PulseData/BaseCalls' not in self.h5handle:
            raise IOError(
                "Unable to locate BaseCalls within hdf5 file (%s)" % h5FN)
        self._BaseCalls = self.h5handle['PulseData/BaseCalls']
        self._holeToBaseIdx = self._buildHoleToIdx(self._BaseCalls['ZMW'])
        self._holeToZMWIdx = self._buildHoleIdx(self._BaseCalls['ZMW'])
        try:
            self._pulseIndex = self._BaseCalls['PulseIndex'][:]
        except ValueError:
            raise IOError(
                "Unable to process pls.h5 file (%s); No bases were called.")
        self._bases = self._BaseCalls['Basecall']
        self._PulseCallCache = {}
        self._softwareDate = None

        if 'PulseData/PulseCalls' in self.h5handle:
            self._PulseCalls = self.h5handle['PulseData/PulseCalls']
            self._holeToPulseIdx = self._buildHoleToIdx(
                self._PulseCalls['ZMW'])

    def __del__(self):
        if self.h5handle:
            self.h5handle.close()

    def _flatten(self, array):
        return array.reshape(np.multiply.reduce(array.shape))

    def _buildHoleIdx(self, zmwData):
        """Builds a mapping of holeNumber:index, where index correctly
        indexes into ZMW-scale tables."""
        holeToIdx = {}
        holeNumbers = np.array(zmwData['HoleNumber'][:])
        holeNumbers = self._flatten(holeNumbers)

        for i in range(len(holeNumbers)):
            holeToIdx[int(holeNumbers[i])] = i

        return holeToIdx

    def _buildHoleToIdx(self, zmwData):
        """Given a ZMW data set (either pulse or base), returns the corresponding hole to index mapping."""
        holeNumbers = np.array(zmwData['HoleNumber'][:])
        holeNumbers = self._flatten(holeNumbers)
        numEvents = np.array(zmwData['NumEvent'][:])
        numEvents = self._flatten(numEvents)

        s = 0
        holeToIdx = {}
        for i in xrange(len(holeNumbers)):
            hn = holeNumbers[i]
            ne = numEvents[i]
            holeToIdx[hn] = (s, ne)
            s += ne
        return holeToIdx

    def _getBaseCallsFromHoleNumber(self, n):
        startIdx, ne = self._holeToBaseIdx[n]
        if ne != 0:
            pulseNumber = self._pulseIndex[startIdx:(startIdx + ne)]
            pulseNumber = self._flatten(pulseNumber)
            basecalls = [chr(c) for c in self._bases[startIdx:(
                startIdx + ne)]]
            #basecalls = list( self._bases[startIdx:(startIdx+ne)] )
        else:
            pulseNumber = []
            basecalls = []
        return pulseNumber, basecalls

    def _getBaseCallsFromXY(self, x, y):
        holeNumber = self._metaData.XY2HoleNumber(x, y)
        return self._getBaseCallsFromHoleNumber(holeNumber)

    def getZMWSeqFromXY(self, x, y):
        dummy, basecalls = self._getBaseCallsFromXY(x, y)
        return "".join(basecalls)

    def hasMetric(self, metric):
        """Returns true iff this pulse file includes the specified
        base-specific metric."""
        if metric in COMPUTED_FIELDS:
            plsDeps = FIELD_DEPENDENCIES[metric][0]
            basDeps = FIELD_DEPENDENCIES[metric][1]
            return all(self.hasMetric(d) for d in plsDeps) or all(self.hasMetric(d) for d in basDeps)
        else:
            return metric in self._BaseCalls or hasattr(self, "_PulseCalls") and metric in self._PulseCalls

    def _getBaseCallAnnotationFromHoleNumber(self, metric, n):
        startIdx, ne = self._holeToBaseIdx[n]
        pulseNumber = self._pulseIndex[startIdx:(startIdx + ne)].transpose()
        if metric not in self._BaseCalls:
            raise MissingMetricError(
                "Unable to locate metric %s in pls.h5 file." % metric)
        values = self._BaseCalls[metric][startIdx:(startIdx + ne)
                                         ] if ne > 0 else np.array([])
        values = self._flatten(values)
        return pulseNumber, values

    def _getPulseFieldForHoleNumber(self, field, n, rpo=True):

        readStartIdx, readNE = self._holeToBaseIdx[n]
        pulseStartIdx, pulseNE = self._holeToPulseIdx[n]

        if field not in self._PulseCallCache:
            self._PulseCallCache[field] = self._PulseCalls[field][:]
        if rpo:
            pulseNumbers = self._pulseIndex[readStartIdx:(
                readStartIdx + readNE)]
            pulseNumbers = self._flatten(pulseNumbers)

            pulseStarts = np.empty(len(pulseNumbers), np.int)
            pulseStarts.fill(pulseStartIdx)
            data = self._PulseCallCache[field][pulseStarts +
                                               pulseNumbers, ...]
        else:
            data = self._PulseCallCache[field][pulseStartIdx:(
                pulseStartIdx + pulseNE)]

        if not np.isscalar(data):
            if min(data.shape) == 1 and len(data.shape) == 2:
                data = self._flatten(data)
        else:
            data = np.array([data])

        return data

    # rpo (Real Pulses Only) should be set to True if you want only
    # pulse data associated with real pulses (called bases).
    def getPulseDataFieldsFromXY(self, fields, x, y, rpo=True):
        holeNumber = self._metaData.XY2HoleNumber(x, y)
        return self.getPulseDataFieldsForHoleNumber(fields, holeNumber, rpo=rpo)

    def _createArrayByChannel(self, array, channelArray):
        if len(array.shape) == 2:
            numberOfChannle = array.shape[1]
        else:  # TODO: fix this where it should be fixed; not here. Hacking this in for AGBT.
            numberOfChannle = array.shape[0]

        # sys.stderr.write(str(array)+"\n\n")
        #sys.stderr.write(str(channelArray) + "\n\n")

        rtnArray = np.zeros(channelArray.shape)
        for i in xrange(numberOfChannle):
            if len(array.shape) == 2:
                rtnArray[channelArray == i] = array[channelArray == i, i]
            else:
                rtnArray[channelArray == i] = [array[i]]

        # sys.stderr.write(str(rtnArray)+"--\n\n")
        return self._flatten(rtnArray)

    def getPulseDataFieldsForHoleNumber(self, fields, hn, rpo=True):

        rtn = {}
        if len(fields) == 0:
            return rtn

        localCache = {}

        # Handle local caching and different locations (base vs pulse)
        def get(field):
            if field not in localCache:
                # if hasattr( self, "_BaseCalls" ) and field in BASECALL_FIELDS
                # and field in self._BaseCalls:
                if hasattr(self, "_BaseCalls") and field in self._BaseCalls:
                    localCache[field] = self._getBaseCallAnnotationFromHoleNumber(
                        field, hn)[1]
                # XXX: JHB - Marting
                elif hasattr(self, "_PulseCalls") and field in self._PulseCalls:
                    localCache[field] = self._getPulseFieldForHoleNumber(
                        field, hn, rpo)
                else:
                    # XXX: This might be risky, but the design is so broken
                    # that it probably doesn't matter.
                    localCache[field] = None
                    # raise MissingMetricError, "Unable to locate metric %s
                    # within the pls file (%s)." % ( field, self.h5FN )
            return localCache[field]

        for field in fields:
            if field not in COMPUTED_FIELDS:
                rtn[field] = get(field)
            else:
                if field == 'StartTime':
                    rtn[field] = get('StartFrame') / self._metaData.FrameRate

                elif field == 'PulseWidth':
                    rtn[field] = get(
                        'WidthInFrames') / self._metaData.FrameRate

                elif field == 'pkmid':
                    rtn[field] = self._createArrayByChannel(
                        get('MidSignal'), get('Channel'))

                elif field == 'pkmax':
                    rtn[field] = self._createArrayByChannel(
                        get('MaxSignal'), get('Channel'))

                elif field == 'pkmean':
                    rtn[field] = self._createArrayByChannel(
                        get('MeanSignal'), get('Channel'))

                elif field == 'IPD' and "PreBaseFrames" in self._BaseCalls:
                    rtn[field] = get(
                        'PreBaseFrames') / self._metaData.FrameRate

                elif field == 'IPD':
                    starts = get('StartFrame')
                    widths = get('WidthInFrames')
                    rtnArray = np.zeros(len(starts), dtype=np.float)
                    rtnArray[1:] = (starts[1:] -
                                    starts[:-1] - widths[:-1]) / self._metaData.FrameRate
                    rtn[field] = self._flatten(rtnArray)

                elif field == 'Light':
                    pkmean = self._createArrayByChannel(
                        get('MeanSignal'), get('Channel'))
                    widths = 1.0 * get(
                        'WidthInFrames') / self._metaData.FrameRate
                    rtn[field] = self._flatten(pkmean * widths)

        return rtn

    def getZMWMetricsForHoleNumber(self, fields, hn):
        """Returns a value from the /BaseCalls/ZMWMetrics/ table."""
        idx = self._holeToZMWIdx[hn]

        # Simple cache for large speedup
        if hasattr(self, "_lastZMWFields") and self._lastZMWFields == fields:
            return dict([(field, self._lastZMWData[field][idx]) for field in fields if field in self._lastZMWData])

        # Cache failure
        self._lastZMWFields = fields
        self._lastZMWData = {}
        for field in fields:
            if ("ZMWMetrics/%s" % field) in self._BaseCalls:
                self._lastZMWData[
                    field] = self._BaseCalls["ZMWMetrics/%s" % field]
            if ("ZMW/%s" % field) in self._BaseCalls:
                self._lastZMWData[field] = self._BaseCalls["ZMW/%s" % field]

        return dict([(field, self._lastZMWData[field][idx]) for field in fields if field in self._lastZMWData])

    def writeFastaFile(self, runcode, filename):
        f = open(filename, 'w')
        movieName = self._metaData.MovieName
        for hn in self._holeToBaseIdx:
            (x, y) = self._metaData.holeNumber2XY(hn)
            dummy, seq = self._getBaseCallsFromHoleNumber(hn)
            if len(seq) == 0:
                continue
            tag = "x%d_y%d_%s_%s" % (x, y, runcode, movieName)
            f.write(">%s|len=%d" % (tag, len(seq)))
            for x in range(0, len(seq), 60):
                f.write("".join(seq[x:x + 60]))
        f.close()


class PLSH5Error(Exception):

    """Generic Error thrown by the high level pls.h5 API"""
    pass


class Movie(object):

    """The top level object of a high-level API for *.pls.h5 files."""

    # This is of course specific to our current setup.
    READ_TO_XY_REGEXP = r'x(\d+)\_y(\d+)\_'

    DRAGON_RABBIT_GLOB = "*Dragon*"

    X_MIN, X_MAX = 1, 34
    Y_MIN, Y_MAX = 2, 94
    PULSE_METRICS = ("StartTime", "PulseWidth", "pkmid", "pkmax",
                     "Channel", "Chi2",
                     "MaxSignal", "MeanSignal", "MidSignal",
                     "StartFrame", "WidthInFrames")
    ZMW_METRICS = ("ReadScore", "HoleStatus", "Productivity", "HQRegionSNR")

    def __init__(self, input_, pulseMetrics=PULSE_METRICS, zmwMetrics=ZMW_METRICS, reportGlob=DRAGON_RABBIT_GLOB):
        """Takes either a path to a .pls.h5 file, or a read name.
        If given a read name, it will glob for the corresponding
        movie file in /mnt/data*. Note that if someone blows away
        the .pls files for your sequence, or replaces them,
        this will break.

        By default, we load in all of the pulse metrics. To save
        time/memory, you can specify a list of the metrics to import.
        MovieHDF5IO.Movie.PULSE_METRICS contains a comprehensive list."""

        self._pulseMetrics = pulseMetrics
        self._zmwMetrics = zmwMetrics

        self._loadMovieFile(input_)
        self._holes = self._mpData._holeToBaseIdx.keys()

    @property
    def channel2base(self):
        return self._mpData.channel2base

    def _loadMovieFile(self, path):
        """Given the path to a movie file, loads it up"""
        self._mpData = MoviePulseData(path)

    def __getitem__(self, key):
        """Returns a ZMW object for the specified read or hole number."""
        if isinstance(key, (int, long)):
            return self._getZMWByHole(key)
        return self._getZMWByName(key)

    def __len__(self):
        """Returns the number of ZMWs in this Movie."""
        return len(self._holes)

    def __iter__(self):
        """Iterates over the reads in this movie. Uses the X_MIN/MAX and Y_MIN/MAX
        constants to determine the X and Y coordinates for all ZMWs."""
        for hole in self._holes:
            yield self._getZMWByHole(hole)

    def _getZMWByHole(self, hole):
        """Returns a ZMW object for the read a tthe given hole number."""
        return ZMW(self._mpData.getPulseDataFieldsForHoleNumber(self._pulseMetrics, hole),
                   self._mpData.getZMWMetricsForHoleNumber(
                       self._zmwMetrics, hole),
                   self, hole)

    def _getZMWByXY(self, x, y):
        """Returns a ZMW object for the read at X,Y."""
        hole = self._mpData._metaData.XY2HoleNumber(x, y)
        return self._getZMWByHole(hole)

    def _getZMWByName(self, readName):
        """Returns a ZMW object for the specified read"""
        readMatch = re.search(Movie.READ_TO_XY_REGEXP, readName)
        if not readMatch:
            raise PLSH5Error("Unable to extract X,Y coordinates from read name (%s)" %
                             readName)
        x, y = int(readMatch.group(1)), int(readMatch.group(2))
        return self._getZMWByXY(x, y)


class ZMW(object):

    """Provides a high level API to pulse data from a single ZMW.
    Use a Movie() to get a ZMW object."""

    def __init__(self, pulseData, zmwData, parent, hole):
        """Should be constructed by accessing someMovie[readName] or someMovie[x,y]."""
        self._pulseData = pulseData
        self._zmwData = zmwData
        self._parent = parent
        self.holeNumber = hole

    @property
    def channel2base(self):
        return self._parent.channel2base

    def __getattr__(self, name):
        """Looks for a ZMW-specific metric matching the attribute name."""
        if name in self._zmwData:
            return self._zmwData[name]
        raise AttributeError("Unable to locate attribute: %s in ZMW" % name)

    def __len__(self):
        """Returns the length of the read, based on a random pulse metric"""
        return int(self._parent._mpData._holeToBaseIdx[self.holeNumber][1])

    def __getitem__(self, key):
        """Can be indexed by an integer to extract all information from a given position
        in the read (ie. read[1] => { "StartTime":200, ... }).
        Or you can index by the pulse metric to return a vector of data for the read
        (ie. read["StartTime"] => [ 200, 230, 560 ]).
        This means that read[1]["StartTime"] == read["StartTime"][1], though the latter
        should be faster, as it follows the structure of the underlying data here."""

        if isinstance(key, basestring):
            if key not in self._parent._pulseMetrics:
                raise PLSH5Error(
                    "Invalid pulse metric used as key to read: (%s)" % key)
            return self._pulseData[key]
        elif isinstance(key, int):
            try:
                return dict([(metric, self._pulseData[metric][key]) for metric in self._pulseData.keys()])
            except IndexError:
                raise IndexError, "Index out of range (%d) for one or more pulse metrics (lengths: %s)" % \
                                  (key, ", ".join(["%s: %d" % (metric, len(self._pulseData[
                                                                           metric])) for metric in self._pulseData.keys()]))
            except KeyError:
                raise KeyError("One or more of these metrics was not found: {k}".format(k=" ".join(self._pulseData.keys())))
        else:
            raise KeyError("Invalid key (%s) specified as read index. Must be either an int or a pulse metric." % str(key))


def getNumSequencingZMWs(h5file):
    """Much faster utility method for simply
    retrieving the number of ZMWs and the number of HoleStatus==0 ZMWs
    Returns (totZMWs, totSequencingZMWs)
    """
    totZMWs = 0
    totSequencing = 0
    h5 = None
    try:
        h5 = h5py.File(h5file, 'r')
        try:
            status = h5['/PulseData/BaseCalls/ZMW/HoleStatus'].value
            totSequencing = sum(status == 0)
        except Exception, e:
            status = h5['/PulseData/BaseCalls/ZMW/HoleNumber'].value
            totSequencing = len(status)
        totZMWs = len(status)
    finally:
        if h5:
            h5.close()
    return totZMWs, totSequencing


def getMovieNameFromPlsH5(plsH5FileName):
    """
    Extract the Movie name from a pls.h5

    :param plsH5FileName: pls.h5 file name
    :return: (str) movie name

    :note MK: I am following the style of the existing API, which needs to be
    rewritten to provide a easier to use interface.
    """
    f = h5py.File(plsH5FileName, 'r')

    groupName = '/ScanData/RunInfo'

    runInfoGroup = f[groupName]

    attrs = {k: v for k, v in runInfoGroup.attrs.iteritems()}

    movieName = None
    if 'MovieName' in attrs:
        # FIXME There are some new Astro movie names which are listed as
        # lists of movie names with only one movie name
        # This may not be a reliable solution!!
        if isinstance(attrs['MovieName'], (np.ndarray, list)):
            movieName = attrs['MovieName'][0]
        elif isinstance(attrs['MovieName'], basestring):
            movieName = attrs['MovieName']
        else:
            raise Exception("Unsupported Movie Name value")
    return movieName


def getSoftwareVersion(h5file):
    """Much faster utility method for
    retrieving the software name and version (pulse2base) which generated
    this file
    """
    DEFAULT_SOURCE = 'PulseToBase'
    groupName = "/PulseData/BaseCalls"
    if h5file.endswith("ccs.h5"):
        DEFAULT_SOURCE = "PulseToCCS"
        groupName = "/PulseData/ConsensusBaseCalls"
    h5 = h5py.File(h5file, 'r')
    version = h5[groupName].attrs['ChangeListID']

    if 'SoftwareVersion' in h5[groupName].attrs:
        source = h5[groupName].attrs['SoftwareVersion']
    else:
        source = DEFAULT_SOURCE

    # this attribute seems to be 0-length these days
    if len(source) == 0:
        source = DEFAULT_SOURCE
    else:
        source = source.replace(' ', '_')
    h5.close()
    return source, version


def test():
    plsFileName = "/mnt/data/vol4/1850238/0008/Reports for SideshowBob/m100316_082405_Geo_p1_b20.pls.h5"
    #     plsFileName = sys.argv[1]
    mmd = MovieMetaData(plsFileName)
    # print mmdF._metaData
    print mmd._findAttrib("ChipArray", mmd._metaData, "", [])
    print mmd.LaserIntensity
    print mmd.XY2HoleNumber(10, 10)
    print mmd.MovieName

    # Test for basecall # >= pulsecall #
    #     pulseData = MoviePulseData(plsFileName)
    # basedata = pulseData.getPulseDataFieldsFromXY(fields=['StartTime'],x=9,y=49,
    # rpo=True)
    # pulsedata = pulseData.getPulseDataFieldsFromXY(fields=['StartTime'],x=9,y=49,
    # rpo=False)
    # print "All pulses = %s, Just basecalls = %s" %
    # (len(pulsedata['StartTime']),len(basedata['StartTime']))

    # print pulseData._getBaseCallsFromXY(10, 20)
    # print pulseData._getBaseCallsFromHoleNumber(4)
    # print pulseData._getPulseFieldForHoleNumber("Channel", 4)
    # print pulseData._getPulseFieldForHoleNumber("Chi2", 4)
    #pulseData.writeFastaFile("test", "test.fa")

    t, s = getNumSequencingZMWs(plsFileName)
    print 'Total Number of ZMWs = %d' % t
    print 'Total Number of Sequencing ZMWs = %d' % s

    s = getSoftwareVersion(plsFileName)
    print 'Software version = %s %s' % s

if __name__ == '__main__':
    test()
