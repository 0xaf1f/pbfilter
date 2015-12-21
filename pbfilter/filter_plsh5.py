#!/usr/bin/env python
"""This script takes in a pls.h5 FOFN (.fofn) plus one or more filter
specifications. A filter specification takes the form: filterName:threshold.
For a list of filters and valid thresholds, use the --availableFilters option.

Also note that trimming by HQRegions is on by default, and can be disabled
using --trim 'False' """
import os
from pprint import pformat
import sys
import optparse
import logging
import random
import argparse

import h5py
import numpy

from pbcore.io import BasH5IO

import pbfilter.io as MovieHDF5IO
from pbfilter.utils import fofn_to_files, validate_fofn, setup_log

log = logging.getLogger(__name__)
__version__ = '2.3'


class PlsRgnFilter(object):

    """Top level class provides a high level API to the user"""

    def __init__(self, plsFofn, filters, rgnFofn=None, plsFilterClass=None):
        """
        Takes a path to a .fofn file of .pls.h5 files, and a
        list of ZMWFilter Objects to filter with
        """
        if rgnFofn is None:
            rgnFofn = plsFofn

        self._plsFNs = fofn_to_files(plsFofn)
        self._rgnFNs = fofn_to_files(rgnFofn)

        self._filters = filters

        self._plsFilterClass = plsFilterClass

        if self._plsFilterClass is None:
            self._plsFilterClass = PlsFilter

    def _writeSummaryHeader(self, summaryOut):
        """Writes out a summary of the filtering to the specified file."""
        headers = ["Movie", "ReadId", "#Bases", "Readlength", "ReadScore"]
        for h in [f.headers for f in self._filters]:
            headers.extend(h)
        headers.append("PassedFilter")
        summaryOut.write(",".join(headers) + "\n")

    def writeFilteredRgnFiles(self, rgnDir, rgnFofn, summaryFN=None, trim=True):
        """Given a location to write the filtered region files and
           an associated FOFN, performs the filtering."""

        if not os.path.isdir(rgnDir):
            os.mkdir(rgnDir)

        summaryOut = open(summaryFN, 'w') if summaryFN is not None else None

        if summaryOut:
            self._writeSummaryHeader(summaryOut)

        self._plsFilters = [self._plsFilterClass(plsFN, rgnFN, self._filters)
                            for plsFN, rgnFN in zip(self._plsFNs, self._rgnFNs)]

        for pf in self._plsFilters:
            pf.writeFilteredRgnFile(rgnDir, trim)

            if summaryOut:
                summaryOut.write("\n".join([",".join(s) for s in pf.summaries]) + '\n')

            pf.clear()  # This can save Gb of memory for large datasets

        if summaryOut:
            summaryOut.close()

        with open(rgnFofn, 'w') as f:
            rgns = [pf.newRgnPath for pf in self._plsFilters]
            log.info("Writing new regions to {f}".format(f=rgnFofn))
            log.debug(rgns)
            f.write("\n".join(rgns) + "\n")


class PlsFilter(object):

    """Filters a single pls/rgn file pair into a single rgn file."""

    def __init__(self, plsFN, rgnFN, filters):
        """Parameters are a plsFN, a rgnFN and a
           list of ZMWFilters to filter with."""
        self._plsFN = plsFN
        self._filters = filters
        self.summaries = []
        self._rgnFN = rgnFN
        self._rbh = None
        self._hqrt = None
        self._newRgnFN = None
        try:
            self._types = list(MovieHDF5IO.regionTypes(rgnFN))
            self._descs = list(MovieHDF5IO.regionDescriptions(rgnFN))
            self._srcs = list(MovieHDF5IO.regionSources(rgnFN))
        except IOError, e:
            raise IOError(e.errno, "Unable to process H5 file: %s\n%s" % (rgnFN, e.strerror))

    def clear(self):
        """Call after all processing is finished and summaries have been extracted.
        Frees up memory."""
        del self.summaries
        del self._baxH5Reader

    @property
    def _rgnsByHole(self):
        """Returns a dictionary mapping hole numbers to lists of regions"""
        if self._rbh is None:
            self._rbh = {}
            for region in MovieHDF5IO.regionIterator(self._rgnFN):
                hn = region["HoleNumber"]
                self._rbh.setdefault(hn, [])
                self._rbh[hn].append(region)
        return self._rbh

    @property
    def _hqRgnType(self):
        """Returns the Type ID for HQRegions. If one does not exist,
        one is created and its ID is returned."""
        if self._hqrt is None:
            if "HQRegion" not in self._types:
                self._types.append("HQRegion")
                self._descs.append("High Quality Region of Sequence")
                self._srcs.append("Filter Module")
            self._hqrt = self._types.index("HQRegion")

        return self._hqrt

    def writeFilteredRgnFile(self, rgnDir, trim):
        """Performs the actual filtering using the stored pls.h5 and filters."""

        try:
            self._baxH5Reader = BasH5IO.BaxH5Reader(self._plsFN)
        except IOError, e:
            raise IOError(e.errno, "Unable to process pls file: %s.\n%s"
                          % (self._plsFN, e.strerror))

        regions = []
        holeNumbers = self._baxH5Reader._mainBasecallsGroup["ZMW/HoleNumber"].value
        for holeNumber in holeNumbers:
            zmw = self._baxH5Reader[holeNumber]
            regions.extend(self._processZMW(zmw, trim))

        regions.sort(key=lambda r: r["HoleNumber"])

        rgnFN = self._plsFN.replace(".pls.h5", ".rgn.h5").replace(".plx.h5", ".rgn.h5").replace(".bas.h5", ".rgn.h5").replace(".bax.h5", ".rgn.h5")
        self._newRgnFN = os.path.join(rgnDir, os.path.basename(rgnFN))

        MovieHDF5IO.writeRegionsTable(regions, self._newRgnFN,
                                      types=self._types,
                                      descriptions=self._descs,
                                      sources=self._srcs)

        log.info("completed writing rgn {f}".format(f=self._newRgnFN))

        plsh5 = h5py.File(self._plsFN, 'r')
        rgnh5 = h5py.File(self._newRgnFN, 'a')

        scanGroupName = '/ScanData/RunInfo'

        rgnh5.create_group('/ScanData')
        rgnRunInfoGroup = rgnh5['/ScanData']

        #log.info("writing '/ScanData/RunInfo to rgn.h5 file {f}".format(f=self._newRgnFN))
        msg = "writing /ScanData/RunInfo to rgn.h5 file '{f}'".format(f=self._newRgnFN)
        log.info(msg)
        plsh5.copy(scanGroupName, rgnRunInfoGroup)

        plsh5.close()
        rgnh5.close()
        log.info("Completed writing region file {r}".format(r=self._newRgnFN))

    def _processZMW(self, zmw, trim):
        """Applies all filters to the specified ZMW, returning the resulting output regions."""
        
        hqRgn = None
        rgns = []
        for regionRow in zmw.regionTable:
            rgn = MovieHDF5IO.PlsRegion(
                HoleNumber=zmw.holeNumber, Start=regionRow['regionStart'],
                End=regionRow['regionEnd'], Score=regionRow['regionScore'],
                TypeIndex=regionRow['regionType'])
            rgns.append(rgn)
            if rgn['TypeIndex'] == BasH5IO.HQ_REGION:
                hqRgn = rgn
                                        
        offsets = zmw.baxH5._offsetsByHole[zmw.holeNumber]
        zmwLength = offsets[1] - offsets[0]
        if hqRgn is None:
            hqRgn = MovieHDF5IO.PlsRegion(HoleNumber=zmw.holeNumber,
                                          TypeIndex=BaseMap.HQ_REGION,
                                          Start=0,
                                          End=zmwLength,
                                          Score=1000)
            if hasattr(zmw, "readScore"):
                hqRgn["Score"] = zmw.readScore
            rgns.append(hqRgn)
        
        if not trim:
            hqRgn["Start"] = 0
            hqRgn = zmwLength
        movie = zmw.baxH5.movieName
        read = zmw.zmwName
        nBases = zmwLength
        readlength = hqRgn["End"] - hqRgn["Start"]
        readScore = hqRgn["Score"] / 1000.0 if hqRgn["Score"] > 1 else hqRgn["Score"]
        summary = [movie, read, "%d" % nBases, "%d" % readlength, "%.4f" % readScore]
        for filtr in self._filters:
            summary.extend(filtr.apply(movie, zmw, hqRgn, rgns, self._types))
        passed = hqRgn["End"] - hqRgn["Start"] > 0
        summary.append("1" if passed else "0")
        self.summaries.append(summary)
        return rgns

    @property
    def rgnPath(self):
        return self._rgnFN

    @property
    def newRgnPath(self):
        return self._newRgnFN


class IZMWFilter(object):

    """Defines the interface to ZMWFilter classes."""

    def __init__(self, headers=(), pulseMetrics=()):
        """
        Defines the column headers for any summary metrics returned by
        apply
        """
        self.headers = headers
        self.pulseMetrics = pulseMetrics

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        """Takes in a path to a pls.h5 file, ZMW object, a Region object
        representing the HQRegion, and a list of all regions for this ZMW.
        Should modify the regions in place and return a list of summary metrics."""
        return []

    def _fail(self, hqRgn):
        """Simple helper function to fail a ZMW"""
        hqRgn["Start"] = 0
        hqRgn["End"] = 0

    def __str__(self):
        h = 'headers=' + ",".join(self.headers) if self.headers else ""
        m = 'metrics=' + ",".join(self.pulseMetrics) if self.pulseMetrics else ""
        return "{k} {h} {m}".format(k=self.__class__.__name__, h=h, m=m)

    def __repr__(self):
        return "<{k}>".format(k=str(self))


class MinReadLengthFilter(IZMWFilter):

    """Filters on a minimum read length cut"""

    def __init__(self, minRL=50):
        super(MinReadLengthFilter, self).__init__()
        self._minRL = float(minRL)

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        trimmedLength = hqRgn["End"] - hqRgn["Start"]
        if trimmedLength < self._minRL:
            self._fail(hqRgn)
        return []


class MaxReadLengthFilter(IZMWFilter):

    """Filters on a maximum read length cut"""

    def __init__(self, maxRL=3000):
        super(MaxReadLengthFilter, self).__init__()
        self._maxRL = float(maxRL)

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        trimmedLength = hqRgn["End"] - hqRgn["Start"]
        if trimmedLength > self._maxRL:
            self._fail(hqRgn)
        return []


class SequencingZMWFilter(IZMWFilter):

    """Filters out ZMWs not marked as sequencing ZMWs"""

    def __init__(self):
        super(SequencingZMWFilter, self).__init__(headers=["SequencingZMW"])

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        sequencingZMW = True
        if zmw.holeNumber not in zmw.baxH5.allSequencingZmws:
            sequencingZMW = False
        #if hasattr(zmw, 'HoleStatus'):
        #    sequencingZMW = zmw.HoleStatus == 0
        if not sequencingZMW:
            self._fail(hqRgn)
        return ["1" if sequencingZMW else "0"]


class ReadScoreFilter(IZMWFilter):

    """Filters on ReadScore"""

    def __init__(self, minReadScore=0.7):
        super(ReadScoreFilter, self).__init__(headers=["Productivity"])
        self._minReadScore = float(minReadScore)

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        """Filters HQRegions on ReadScore"""

        if hqRgn["Score"] > 1.0:
            score = hqRgn["Score"] / 1000.0
        else:
            score = hqRgn["Score"]

        if score < self._minReadScore:
            self._fail(hqRgn)

        return ["%d" % (zmw.productivity if hasattr(zmw, "productivity") else -1)]


class SNRFilter(IZMWFilter):

    def __init__(self, minSNR=0.0):
        super(SNRFilter, self).__init__(headers=["SNR"])
        self._minSNR = float(minSNR)

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        #snr = min(zmw.HQRegionSNR)
        snr = min(zmw.zmwMetric("HQRegionSNR"))
        if snr < self._minSNR:
            self._fail(hqRgn)
        return ["%.4f" % snr]


class InsertAdapterFilter(IZMWFilter):

    """Filters Inserts and Adapters by the HQRegion."""

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        toRemove = []
        trimTypes = [BasH5IO.ADAPTER_REGION, BasH5IO.INSERT_REGION]
        #trimTypes = [rgnTypes.index(t) for t in ["Insert", "Adapter"]
        #             if t in rgnTypes]
        if len(trimTypes) > 0:
            for rgn in rgns:
                if rgn["TypeIndex"] in trimTypes:
                    rgn["Start"] = max(rgn["Start"], hqRgn["Start"])
                    rgn["End"] = min(rgn["End"], hqRgn["End"])
                if (rgn["End"] - rgn["Start"] <= 0 and
                        rgn["TypeIndex"] != BasH5IO.HQ_REGION):
                    toRemove.append(rgn)
        for rgnToRemove in toRemove:
            rgns.remove(rgnToRemove)
        return []


class ReadIdFilter(IZMWFilter):

    """Accepts all reads that occur in the specified list of 'white' reads."""

    def __init__(self, whitelistFN):
        super(ReadIdFilter, self).__init__(headers=["Whitelisted"])
        # set of tuples {(movie_name, hole_number), ...}
        self._whitelist = set()
        self._whitelist_file = whitelistFN
        try:
            with open(self._whitelist_file, 'r') as f:
                for line in f:
                    # default error message
                    emsg = "Invalid whitelist format for '{s}' in file {f}.".format(s=line.strip(), f=self._whitelist_file)
                    s = line.strip().split('/')
                    if len(s) != 2:
                        log.error(emsg)
                        raise ValueError(emsg)
                    else:
                        # this should probably be validated too
                        movie_name = s[0]
                        try:
                            hole_number = int(s[1])
                        except TypeError:
                            m = "Invalid hole number type {t}.".format(t=s[1]) + emsg
                            log.error(m)
                            raise TypeError(m)

                        self._whitelist.add((movie_name, hole_number))

        except IOError as e:
            e.strerror = "Unable to find in ReadId Whitelist (%s)" % whitelistFN
            log.error(e)
            raise e

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        white = 1
        if (movie, zmw.holeNumber) not in self._whitelist:
            self._fail(hqRgn)
            white = 0
        return ["%d" % white]


class SubsamplingFilter(IZMWFilter):

    """Randomly accepts reads at a rate acceptRate <= 1.0"""

    def __init__(self, acceptRate=1.0):
        super(SubsamplingFilter, self).__init__()
        self._acceptRate = float(acceptRate)

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        p = random.random()
        if p >= self._acceptRate:
            self._fail(hqRgn)
        return []


class ISubreadFilter(IZMWFilter):

    """Base class for Subread filters"""

    def __init__(self):
        super(ISubreadFilter, self).__init__()

    def _regionFailsFilter(self, rgn):
        """Does this region pass the filter"""
        return False

    def apply(self, movie, zmw, hqRgn, rgns, rgnTypes):
        insertIdx = BasH5IO.INSERT_REGION
        toRemove = []
        for rgn in rgns:
            if rgn["TypeIndex"] == insertIdx and self._regionFailsFilter(rgn):
                toRemove.append(rgn)
        for rgnToRemove in toRemove:
            rgns.remove(rgnToRemove)
        return []


class MinSubreadLengthFilter(ISubreadFilter):

    """Removes all insert regions less than min subread length."""

    def __init__(self, minSRL=50):
        super(MinSubreadLengthFilter, self).__init__()
        self._minSRL = float(minSRL)

    def _regionFailsFilter(self, rgn):
        return (rgn["End"] - rgn["Start"]) < self._minSRL


class MaxSubreadLengthFilter(ISubreadFilter):

    """Removes all insert regions more than max subread length."""

    def __init__(self, maxSRL=3000):
        super(MaxSubreadLengthFilter, self).__init__()
        self._maxSRL = float(maxSRL)

    def _regionFailsFilter(self, rgn):
        return (rgn["End"] - rgn["Start"]) > self._maxSRL


SUPPORTED_FILTERS = {"MinReadScore": ("Accuracy Prediction from Primary",
                                      "0.0-1.0 (e.g. 0.75)",
                                      ReadScoreFilter),
                     "MinRL": ("Minimum # of Bases in the ZMW",
                               "0.0-Inf (e.g. 50)",
                               MinReadLengthFilter),
                     "MaxRL": ("Maximum # of Bases in the ZMW",
                               "1.0-Inf (e.g. 10000)",
                               MaxReadLengthFilter),
                     "ReadWhitelist": ("Newline separated list of ReadIds",
                                       "file_path",
                                       ReadIdFilter),
                     "Subsampling": ("Randomly subsample reads",
                                     "0.0-1.0",
                                     SubsamplingFilter),
                     "MinSRL": ("Minimum # of Bases in a Subread",
                                "0.0-Inf (e.g. 50)",
                                MinSubreadLengthFilter),
                     "MaxSRL": ("Maximum # of Bases in a Subread",
                                "0.0-Inf (e.g. 10000)",
                                MaxSubreadLengthFilter)
                     }


def _print_available_filters(supported_filters):
    """Prints information on available filters and their thresholds."""
    widths = (20, 40, 20)
    data = [("Filter", "Description", "Threshold Values"),
            ("------", "-----------", "----------------")]
    # this is stupid
    for f, (d, t, c) in supported_filters.items():
        data.append((f, d, t))

    print

    for row in data:
        i = 1
        nextline = "\n"
        for col, width in zip(row, widths):
            print col[:width] + " " * max(0, width - len(col)),
            if not i == 2:
                i += 1
                continue
            mycol = col[width:]
            mybgn = width + 1
            while len(mycol) > 1:
                nextline += " " * 21
                nextline += mycol[:width]
                nextline += " " * (width - len(mycol))
                nextline += "\n"
                mycol = mycol[width:]
                mybgn += width
            i += 1
        print nextline,
    print
    return 0


def get_parser():
    parser = argparse.ArgumentParser(version=__version__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_fofn', help="Path to bas.fofn",
                        type=validate_fofn)

    parser.add_argument("--outputSummary", default=None,
                        help="File to write filtering information to.")

    parser.add_argument("--outputDir", type=str, default=os.getcwd(),
                        help="Directory to write regions table results (rgn.h5s).")

    parser.add_argument("--outputFofn", type=str, default="filtered_regions.fofn",
                        help="A FOFN file pointing to the regions files.")

    parser.add_argument("--debug", action="store_true",
                        help="Outputs a log to stderr with helpful debug info.")

    parser.add_argument("--logFile", default=None,
                        help="Set a log file for logging output. Defaults to stderr.")

    # this needs to be fixed.
    parser.add_argument("--trim", default='True',
                        help="Set to 'True' to trim down to HQRegion. Default is True.")

    # this should validate at the argparse level
    parser.add_argument("--filter", default="",
                        help="Define the set of filters to be run. "
                        "Ex. --filter='minRL=50,minReadScore=0.7'")

    parser.add_argument("--availableFilters", action="store_true",
                        help="Output descriptions of available filters and exit.")

    return parser


def _to_filters(filter_str):
    """Convert filter str to Filter instances"""
    filters = []
    for spec in filter_str.split(","):
        if len(spec) > 0:
            filterName, threshold = spec.split("=")
            # Grab out the class and initialize it with the threshold.
            filters.append(SUPPORTED_FILTERS[filterName][2](threshold))
    return filters


def run(pls_fofn, output_fofn, output_dir, filters,
        output_summary=None, trim=True, debug=False):
    """Main function to run filters

    Convert these fofns, to use lists of files
    """
    assert(all(isinstance(f, IZMWFilter) for f in filters))

    bas_files = fofn_to_files(pls_fofn)
    log.info("Processing {n} files in {f}:".format(n=len(bas_files), f=pls_fofn))
    log.debug(pformat(bas_files))
    log.info("{n} filters:".format(f=filters, n=len(filters)))
    log.debug(pformat(filters))

    plsRgnFilter = PlsRgnFilter(pls_fofn, filters)

    try:
        plsRgnFilter.writeFilteredRgnFiles(output_dir, output_fofn,
                                           summaryFN=output_summary,
                                           trim=trim)
    except IOError as e:
        log.error("I/O Error in accessing pls.h5 files. (%s)" % pls_fofn)
        log.error(e)
        sys.stderr.write(str(e) + "\n")
        return 1

    return 0


def main():
    p = get_parser()

    args = p.parse_args()
    debug = args.debug
    input_fofn = args.input_fofn
    output_dir = args.outputDir
    output_summary_csv = args.outputSummary
    output_fofn = args.outputFofn
    log_file = args.logFile
    trim = args.trim
    filter_str = args.filter
    available_filters = args.availableFilters

    if available_filters:
        return _print_available_filters(SUPPORTED_FILTERS)

    filters = []
    # These are run by default.
    filters.append(SequencingZMWFilter())
    filters.append(InsertAdapterFilter())

    filters.extend(_to_filters(filter_str))

    if debug:
        setup_log(log, level=logging.DEBUG, file_name=log_file)
    else:
        log.addHandler(logging.NullHandler())

    rcode = run(input_fofn, output_fofn, output_dir, filters,
                output_summary=output_summary_csv, trim=trim, debug=debug)

    log.info("Exiting {f} v{v} with return code {r}.".format(r=rcode, f=os.path.basename(__file__), v=__version__))
    return rcode
