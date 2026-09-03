import ast
import gzip
import re
import h5py
import numpy
from famdb_globals import (
    SOUNDEX_LOOKUP,
    GROUP_FAMILIES,
    dfam_acc_pat,
)
from famdb_helper_classes import Family


def accession_bin(acc):
    """Maps an accession (Dfam or otherwise) into apropriate bins (groups) in HDF5"""
    dfam_match = dfam_acc_pat.match(acc)
    if dfam_match:
        path = (
            GROUP_FAMILIES
            + "/"
            + dfam_match.group(1)
            + "/"
            + dfam_match.group(2)
            + "/"
            + dfam_match.group(3)
            + "/"
            + dfam_match.group(4)
        )
    else:
        path = GROUP_FAMILIES + "/Aux/" + acc[0:2].lower()
    return path


def get_family(entry):
    """Builds a Family from db data"""
    if not entry:
        return None

    family = Family()

    # Read the family attributes and data
    for k in entry.attrs:
        value = entry.attrs[k]
        if k == "model":
            # Normalise whatever storage format was used for the model so that
            # callers always receive a plain decoded string.
            #
            # Historical storage variants (all produced by intermediate export
            # code before the current sibling-dataset format was stabilised):
            #
            #   (a) Raw gzip bytes stored as HDF5 bytes attribute.
            #       h5py returns bytes / numpy.bytes_.
            #
            #   (b) Python repr of gzip bytes stored as a string attribute,
            #       e.g. the literal text "b'\x1f\x8b...'" - produced when
            #       str(compressed_bytes) was accidentally called before
            #       storing.  h5py returns a str starting with "b'".
            if isinstance(value, (bytes, numpy.bytes_)):
                value = gzip.decompress(bytes(value)).decode()
            elif isinstance(value, str) and value.startswith("b'"):
                try:
                    raw = ast.literal_eval(value)   # str -> bytes
                    value = gzip.decompress(raw).decode()
                except Exception:
                    pass  # leave as-is; to_dfam_hmm will handle or skip
        setattr(family, k, value)

    # Preferred storage (current format): gzip-compressed sibling dataset.
    # This overrides anything that may have been read from attrs above.
    model_key = entry.name.split("/")[-1] + ".model"
    if model_key in entry.parent:
        compressed = entry.parent[model_key][()].tobytes()
        decoded = gzip.decompress(compressed).decode()
        # Handle the case where str(gzip_bytes) was accidentally stored instead
        # of the raw bytes, producing gzip(repr(gzip(HMMER3))).  One layer of
        # gzip is already stripped above; if the result looks like a bytes repr
        # we use ast.literal_eval to recover the actual bytes and decompress once
        # more to reach the original HMMER3 text.
        if decoded.startswith("b'") or decoded.startswith('b"'):
            try:
                inner = ast.literal_eval(decoded)  # repr-string -> bytes
                decoded = gzip.decompress(inner).decode()
            except Exception:
                pass  # leave as-is; to_dfam_hmm will handle or skip
        family.model = decoded

    return family


def families_iterator(g, prefix=""):
    """Generator that returns all family accession keys in a group.

    Skips the companion '.model' datasets that are stored alongside each
    family dataset in the current file format.
    """
    for key, item in g.items():
        path = f"{prefix}/{key}"
        if isinstance(item, h5py.Dataset):  # test for dataset
            if not key.endswith(".model"):   # skip compressed-model companions
                yield (key)
        elif isinstance(item, h5py.Group):  # test for group (go down)
            yield from families_iterator(item, path)


# Filter methods --------------------------------------------------------------------------
def filter_name(family, name):
    """Returns True if the family's name begins with 'name'."""

    if family.attrs.get("name"):
        if family.attrs["name"].lower().startswith(name):
            return True

    return False


def filter_search_stages(family, stages):
    """Returns True if the family belongs to a search stage in 'stages'."""
    if family.attrs.get("search_stages"):
        sstages = (ss.strip() for ss in family.attrs["search_stages"].split(","))
        for family_ss in sstages:
            if family_ss in stages:
                return True

    return False


# RMH: 6/27/25
def filter_defined_search_stages(family):
    """Returns True if the family has search stages defined."""
    if family.attrs.get("search_stages"):
        return False

    return True


def filter_repeat_type(family, rtype):
    """
    Returns True if the family's RepeatMasker Type plus SubType
    (e.g. "DNA/CMC-EnSpm") starts with 'rtype'.
    """
    if family.attrs.get("repeat_type"):
        full_type = family.attrs["repeat_type"]
        if family.attrs.get("repeat_subtype"):
            full_type = full_type + "/" + family.attrs["repeat_subtype"]

        if full_type.lower().startswith(rtype):
            return True

    return False


def filter_curated(accession, curated):
    """
    Returns True if the family's curatedness is the same as 'curated'. In
    other words, 'curated=True' includes only curated familes and
    'curated=False' includes only uncurated families.

    Families are currently assumed to be curated unless their name is of the
    form DR<9 digits>.

    TODO: perhaps this should be a dedicated 'curated' boolean field on Family
    """

    is_curated = not (
        accession.startswith("DR")
        and len(accession) == 11
        and all((c >= "0" and c <= "9" for c in accession[2:]))
    )

    return is_curated == curated


def soundex(word):
    """
    Converts 'word' according to American Soundex[1].

    This is used for "sounds like" types of searches.

    [1]: https://en.wikipedia.org/wiki/Soundex#American_Soundex
    """

    codes = [SOUNDEX_LOOKUP[ch] for ch in word.upper() if ch in SOUNDEX_LOOKUP]

    # Start at the second code
    i = 1

    # Drop identical sounds and H and W
    while i < len(codes):
        code = codes[i]
        prev = codes[i - 1]

        if code is None:
            # Drop H and W
            del codes[i]
        elif code == prev:
            # Drop adjacent identical sounds
            del codes[i]
        else:
            i += 1

    # Keep the first letter
    coding = word[0]

    # Keep codes, except for the first or vowels
    codes_rest = filter(lambda c: c > 0, codes[1:])

    # Append stringified remaining numbers
    for code in codes_rest:
        coding += str(code)

    # Pad to 3 digits
    while len(coding) < 4:
        coding += "0"

    # Truncate to 3 digits
    return coding[:4]


def sounds_like(first, second):
    """
    Returns true if the string 'first' "sounds like" 'second'.

    The comparison is currently implemented by running both strings through the
    soundex algorithm and checking if the soundex values are equal.
    """
    soundex_first = soundex(first)
    soundex_second = soundex(second)

    return soundex_first == soundex_second


def sanitize_name(name):
    """
    Returns the "sanitized" version of the given 'name'.
    This must be kept in sync with Dfam's algorithm.
    """
    name = re.sub(r"[\s\,\_]+", "_", name)
    name = re.sub(r"[\(\)\<\>\']+", "", name)
    return name


def is_fasta(infile):
    fasta_el = {"header": None, "body": None}
    with open(infile, "r") as file:
        for line in file.readlines():

            if line.startswith(">") and fasta_el["header"] is not None:
                fasta_el["header"] = line
            elif not line.startswith(">") and fasta_el["body"] is not None:
                fasta_el["body"] = line

            if fasta_el["header"] is not None and fasta_el["body"] is not None:
                fasta_el["header"] = None
                fasta_el["body"] = None
    return fasta_el["header"] is None and fasta_el["body"] is None