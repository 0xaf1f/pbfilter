import sys
import logging
import tempfile
import string
from subprocess import CalledProcessError

from pbcore.util.Process import backticks
from pbcore.io.FastaIO import FastaReader, FastaRecord, FastaWriter

# External Exe
from pbfilter.utils import log_timing

_PLS_TO_FASTA_EXE = 'pls2fasta'
_FASTA_REV_COMP_EXE = 'fastarevcomp'
_SPD_MATCHER_EXE = 'sdpMatcher'

_COMPLEMENT_TRANSFORM = string.maketrans('acgtrymkbdhvACGTRYMKBDHV',
                                         'tgcayrkmvhdbTGCAYRKMVHDB')

log = logging.getLogger(__name__)


def _sequence_to_reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT_TRANSFORM)[::-1]


def fasta_to_reverse_complement(fasta_input, fasta_output):
    """Write the reverse complement of a fasta file to a fasta file"""

    with FastaReader(fasta_input) as r:
        with FastaWriter(fasta_output) as w:
            for record in r:
                rc_sequence = _sequence_to_reverse_complement(record.sequence)
                reverse_record = FastaRecord(record.name, rc_sequence)
                w.writeRecord(reverse_record)

    # follow the similar subprocess model
    return 0


def _run_backticks(cmd):
    log.info("Running cmd:")
    log.info(cmd)

    stdout, rcode, stderr = backticks(cmd)

    if rcode != 0:
        msg = "Non-zero ({r}) for cmd {c}".format(c=cmd, r=rcode)
        log.error(msg)
        log.error(stdout)
        log.error(stderr)
        sys.stderr.write(msg + "\n")
        sys.stderr.write(stderr + "\n")
        raise CalledProcessError(rcode, cmd)

    log.debug("Successfully completed command {c}.".format(c=cmd))
    return rcode


@log_timing
def run_pls_to_fasta(pls_file_name, rgn_file_name, fasta_output):
    cmd = _to_pls_to_fasta_cmd(pls_file_name, rgn_file_name, fasta_output)

    return _run_backticks(cmd)


def _to_pls_to_fasta_cmd(pls_file_name, rgn_file_name, fasta_output):
    cmd = "{e} '{p}' '{f}' -regionTable '{r}'".format(e=_PLS_TO_FASTA_EXE,
                                                      p=pls_file_name,
                                                      f=fasta_output,
                                                      r=rgn_file_name)
    return cmd


@log_timing
def _parse_scores(scores_output):
    """Returns a dict of the summary scores {id:score}"""

    with open(scores_output, 'r') as f:
        summary_lines = [l.strip().split(",") for l in f]

    def fields_to_score(f):
        return "/".join(f[0].split("/")[0:2]), int(f[-1])

    summary = {}
    for l in summary_lines[1:]:
        id_, score = fields_to_score(l)
        summary[id_] = min(score, summary.get(id_, 0))

    return summary


@log_timing
def run_sdp_matcher(fasta_input, k, fasta_output, scores_output):

    cmd = _to_sdp_matcher_cmd(fasta_input, k, fasta_output, scores_output)

    return _run_backticks(cmd)


def _to_sdp_matcher_cmd(fasta_file, fasta_revcomp_file, k, scores_output):

    cmd = "{e} {f} {r} {k} -local > {s}".format(e=_SPD_MATCHER_EXE,
                                                f=fasta_file,
                                                r=fasta_revcomp_file,
                                                k=k,
                                                s=scores_output)
    return cmd


def _to_fasta_rev_comp(fasta_input, fasta_output):
    # TODO : Rewrite this in pure python, or cython.
    return "{e} {i} > {o}".format(e=_FASTA_REV_COMP_EXE, i=fasta_input,
                                  o=fasta_output)


@log_timing
def run_fasta_rev_comp_with_exe(fasta_input, fasta_output):
    """A pure python implementation of the reverse complement"""
    return fasta_to_reverse_complement(fasta_input, fasta_output)


@log_timing
def run_fasta_rev_comp(fasta_input, fasta_output):
    cmd = _to_fasta_rev_comp(fasta_input, fasta_output)
    return _run_backticks(cmd)


def _tmp_file(suffix, delete=False):
    tf = tempfile.NamedTemporaryFile(suffix=suffix, delete=delete)
    tf.close()
    log.info("Created tmp file {t}".format(t=tf.name))
    return tf.name


@log_timing
def run_analysis(pls_file, rgn_file, k=10):
    """Return scores dict

    :raises: CalledProcessError if non-zero exit code occurs
    :raises: IOError, TypeError, IndexError if scores can not be parsed
    """

    tmp_fasta = _tmp_file('.fasta')
    rcode = run_pls_to_fasta(pls_file, rgn_file, tmp_fasta)

    tmp_rev_fasta = _tmp_file('reverse_comp.fasta')
    rcode = run_fasta_rev_comp(tmp_fasta, tmp_rev_fasta)

    # is this a CSV file?
    scores_csv = _tmp_file(".scores.csv")
    rcode = run_sdp_matcher(tmp_fasta, tmp_rev_fasta, k, scores_csv)

    scores = _parse_scores(scores_csv)

    return scores