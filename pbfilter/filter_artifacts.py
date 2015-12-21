"""This collection of classes allows for extensible filtering of pls.h5 files
using Regions."""
import sys
import os
import logging
import tempfile
import shutil
import argparse
import time

from pbcore.util.Process import backticks
from pbcore.io import BasH5IO

import pbfilter.io as MovieHDF5IO
from pbfilter.filter_plsh5 import PlsRgnFilter, PlsFilter, IZMWFilter
from pbfilter.utils import setup_log, validate_fofn, log_timing
from pbfilter.external_tools import run_analysis

log = logging.getLogger(__name__)

__version__ = '2.2'


class BaseArtifactFilter(PlsFilter):
    """Base class for filtering a single pls/rgn file tuple into a single
    rgn file."""

    def __init__(self, plsFN, rgnFN, filters):
        PlsFilter.__init__(self, plsFN, rgnFN, filters)
        self._art = None
        self.summaries = []
        self._ranAnalysis = False

    @property
    def _artifactRgnType(self):
        """Returns the Type ID for ArtifactRegions. If one does not exist,
        one is created and its ID is returned."""
        if self._art is None:
            if "ArtifactRegion" not in self._types:
                self._types.append("ArtifactRegion")
                self._descs.append("Region of sequence thought to contain a sequencing artifact")
                self._srcs.append("Artifact Filter Module")
            self._art = self._types.index("ArtifactRegion")
        return self._art

    def _runAnalysis(self):
        k = self._filters[0].k
        scores_dct = run_analysis(self._plsFN, self._rgnFN, k=k)
        self.summaryDict = scores_dct
        self._ranAnalysis = True

    def _processZMW(self, zmw, trim):
        """
        Applies this filter to the specified ZMW, returning the resulting
        output regions.

        """

        if not self._ranAnalysis:
            self._runAnalysis()

        rgns = []
        for regionRow in zmw.regionTable:
            rgn = MovieHDF5IO.PlsRegion(
                HoleNumber=zmw.holeNumber, Start=regionRow['regionStart'],
                End=regionRow['regionEnd'], Score=regionRow['regionScore'],
                TypeIndex=regionRow['regionType'])
            rgns.append(rgn)
        hqRgns = filter(lambda rgn: rgn["TypeIndex"] == BasH5IO.HQ_REGION, rgns)
        hqRgn = hqRgns[0] if len(hqRgns) else None

        offsets = zmw.baxH5._offsetsByHole[zmw.holeNumber]
        zmwLength = offsets[1] - offsets[0]
        movie = zmw.baxH5.movieName
        read = zmw.zmwName

        score = "na"
        threshold = self._filters[0].threshold

        if read in self.summaryDict:
            score = "%d" % self.summaryDict[read]
            artifactRegion = MovieHDF5IO.PlsRegion(HoleNumber=zmw.holeNumber,
                                                   TypeIndex=self._artifactRgnType,
                                                   Start=0,
                                                   End=zmwLength,
                                                   Score=int(score))

            if int(artifactRegion["Score"]) <= int(threshold) and trim and hqRgn:
                hqRgn["Start"] = 0
                hqRgn["End"] = 0

            rgns.append(artifactRegion)

        passed = "1" if (hqRgn["End"] - hqRgn["Start"] > 0) else "0"
        summary = [movie, read, score, passed]

        self.summaries.append(summary)

        return rgns


class ArtifactPlsRgnFilter(PlsRgnFilter):

    """Top level class for filtering library artifacts from pls/rgn files."""

    def __init__(self, plsFofn, filters, rgnFofn=None, plsFilterClass=BaseArtifactFilter):
        PlsRgnFilter.__init__(self, plsFofn, filters, rgnFofn, plsFilterClass)

    def _writeSummaryHeader(self, summaryOut):
        headers = ["Movie", "ReadId", "ArtifactScore", "PassedFilter"]
        summaryOut.write(",".join(headers) + "\n")


class ArtifactFilter(IZMWFilter):

    """Just a storage class for the params for filtering artifacts"""

    def __init__(self, threshold=-sys.maxint, k=10):
        IZMWFilter.__init__(self, headers=(), pulseMetrics=())
        self.threshold = threshold
        self.k = k

SUPPORTED_FILTERS = {
    "ForwardReverseScore":
    ("Score for self vs self-rc BLASR SDP alignment",
     "-inf - 0.0",
     ArtifactFilter),
}


def _validate_filters(filter_str):
    """Validate, but no parse filter_str"""
    keyword, threshold = _parse_filter(filter_str)
    return filter_str


def _parse_filter(filter_str):
    algKeyword, threshold = filter_str.split("=")
    if algKeyword not in SUPPORTED_FILTERS:
        raise ValueError("Supported filter algorithms are: %s" % ",".join(SUPPORTED_FILTERS.keys()))

    # raise TypeError if the threshold isn't correct
    int_threshold = int(threshold)
    return algKeyword, int_threshold


def get_parser():
    parser = argparse.ArgumentParser(version=__version__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('pls_fofn', help="input pls FOFN", type=validate_fofn)
    parser.add_argument('rgn_fofn', help="Region FOFN.", type=validate_fofn)

    parser.add_argument("--outputSummary", default=None,
                        help="Optional name of a CSV file to write filtering information to. If it already exists, will append fields to this file.")

    parser.add_argument("--outputDir", default=os.getcwd(), type=str,
                        help="Specify a directory to write results to in the format of regions tables (rgn.h5s).")

    parser.add_argument("--outputFofn", default="filtered_regions.fofn", type=str,
                        help="Specify the location of a FOFN file pointing to the regions files.")

    parser.add_argument("--debug", action="store_true",
                        help="Outputs a log to stderr with helpful debug info.")

    parser.add_argument("--logFile", default=None,
                        help="Set a log file for logging output. Defaults to stderr.")

    parser.add_argument("--filter", type=_validate_filters, default="ForwardReverseScore=-1000",
                        help="The specific artifact filtering algorithm to use (mostly for testing).")

    parser.add_argument("--trim", default="True", type=str,
                        help="Trim hiqh-quality regions based on results of filter.")

    return parser


@log_timing
def _mergeSummaries(temp_output_csv, output_csv, outputDir):

    log.info("Merging CSV summary from {t} to {o}".format(t=temp_output_csv, o=output_csv))

    with open(output_csv, 'r') as f:
        firstLine = f.readline()
        header = firstLine.strip().split(",")

    if len(header) < 3 or header[0] != "Movie" or header[1] != "ReadId" or header[-1] != "PassedFilter":
        msg = "Attempt to merge with an unexpected file {s} header. Got header '{h}', expected Movie, ReadId and PassedFilter to be CSV".format(s=temp_output_csv, h=header)
        log.error(msg)
        sys.stderr.write(msg + "\n")
        raise ValueError(msg)

    outputSummaryFields = len(header)
    tempCsv = tempfile.NamedTemporaryFile(dir=outputDir).name

    # TODO: Rewrite this in pure python.
    cmd = "paste %s %s -d, | cut -f1-%i,%i- -d, > %s; mv %s %s;" % (output_csv, temp_output_csv, outputSummaryFields - 1, outputSummaryFields + 3, tempCsv, tempCsv, output_csv)

    log.info(cmd)
    stdout, rcode, stderr = backticks(cmd)

    if rcode != 0:
        log.fatal("Problem running %s - %s" % (cmd, stderr))
        raise SystemExit

    log.info("completed merging CSV summaries")


def run(pls_fofn, rgn_fofn, output_fofn, output_dir, algKeyword, threshold, summary_csv=None):
    """Actually performs the filtering of artifacts."""

    log.info("Starting run.")
    started_at = time.time()

    description, info, filter_klass = SUPPORTED_FILTERS[algKeyword]
    myFilter = filter_klass(threshold)
    artifactFilter = ArtifactPlsRgnFilter(pls_fofn, [myFilter], rgnFofn=rgn_fofn)

    # Used for merging exiting filtering done in filter_plsh5.py. This should
    # be refactored to have the filtering done at the same step. The filtering
    # shouldn't be done as two distinct operations.
    temp_csv_summary = tempfile.NamedTemporaryFile(dir=output_dir, delete=False)
    temp_csv_summary.close()
    temp_csv_summary_name = temp_csv_summary.name

    try:
        artifactFilter.writeFilteredRgnFiles(output_dir, output_fofn,
                                             summaryFN=temp_csv_summary_name)

    except IOError as e:
        log.error("I/O Error in accessing pls.h5 files (%s) or rgn.h5 files (%s)" % (pls_fofn, rgn_fofn))
        log.error(e)
        return 1

    if summary_csv is not None:
        if os.path.exists(summary_csv):
            _mergeSummaries(temp_csv_summary_name, summary_csv, output_dir)
            os.remove(temp_csv_summary_name)
        else:
            shutil.move(temp_csv_summary_name, summary_csv)

    run_time = time.time() - started_at
    log.info("Successfully completed run in {s:.2f} sec.".format(s=run_time))
    return 0


def main():
    """Main point of entry"""
    p = get_parser()
    args = p.parse_args()

    pls_fofn = args.pls_fofn
    rgn_fofn = args.rgn_fofn
    output_fofn = args.outputFofn
    summary_csv = args.outputSummary
    output_dir = args.outputDir
    debug = args.debug
    filter_str = args.filter

    algo_keyword, threshold = _parse_filter(filter_str)

    if debug:
        setup_log(log, level=logging.DEBUG)
        log.info("Starting {f} v{v}".format(f=os.path.basename(__file__),
                                            v=__version__))
        log.info(args)
    else:
        log.addHandler(logging.NullHandler())

    rcode = run(pls_fofn, rgn_fofn, output_fofn, output_dir, algo_keyword,
                threshold, summary_csv=summary_csv)

    log.info("Exiting {f} v{v}".format(f=os.path.basename(__file__), v=__version__))
    return rcode
