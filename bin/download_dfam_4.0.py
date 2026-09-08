#!/usr/bin/env python3
# 
# adapted from https://github.com/Dfam-consortium/FamDB
"""
download_dfam.py - Interactive downloader for FamDB component files from Dfam.

Fetches the available release from the Dfam server, lets you choose which
components and partitions to download, validates MD5 checksums, and
decompresses the files into Libraries/famdb/ (or a user-specified directory).

Recoverable: already-decompressed files are skipped; partially-downloaded
.gz files are re-verified before re-downloading.

Usage: download_dfam.py [-h] [-o OUTPUT_DIR] [-u URL] [--dry-run]
"""

import argparse
import subprocess
import configparser
import gzip
import hashlib
import os
import re
import shutil
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from html.parser import HTMLParser

DEFAULT_URL = "https://www.dfam.org/releases/current/families/FamDB/"
#DEFAULT_URL = "https://www.dfam.org/releases/Dfam_4.0/families/FamDB/"
# This script lives in <install_dir>/utils/, so go up one level.
_INSTALL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _detect_existing_prefix(out_dir):
    """Return the release prefix found in out_dir's .h5 files, or None."""
    if not os.path.isdir(out_dir):
        return None
    for fname in sorted(os.listdir(out_dir)):
        if ".h5" in fname:
            m = re.match(r"^(dfam\w+?)\.", fname)
            if m:
                return m.group(1)
    return None


# Component display order and friendly names
COMPONENT_ORDER = ["root", "curated.consensus", "curated.hmm", "uncurated.consensus", "uncurated.hmm"]
COMPONENT_LABELS = {
    "root":                 "Root (taxonomy index)",
    "curated.consensus":    "Curated consensus sequences",
    "curated.hmm":          "Curated profile HMMs",
    "uncurated.consensus":  "Uncurated consensus sequences",
    "uncurated.hmm":        "Uncurated profile HMMs",
}

FAMDB_OUTPUT_KEYS_TO_COMPONENTS = {
    'Curated Consensus': "curated.consensus",
    'Curated HMMs': "curated.hmm",
    'Uncurated Consensus': "uncurated.consensus",
    'Uncurated HMMs': "uncurated.hmm"
}


class _LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for attr, val in attrs:
                if attr == "href" and val and not val.startswith(("?", "/", "http")):
                    self.links.append(val.rstrip("/"))

    def error(self, message):
        pass


def _fetch_html(url):
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return resp.read().decode("utf-8")
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach {url}: {e}") from e


def fetch_release_tag(release_root_url):
    """
    Look for a RELEASE_DFAM_X_Y tag file in the release root directory listing.
    Returns (tag, major, minor) or raises RuntimeError.
    """
    html = _fetch_html(release_root_url)

    parser = _LinkParser()
    parser.feed(html)

    release_pat = re.compile(r"^RELEASE_DFAM_(\d+)_(\d+)$")
    for link in parser.links:
        m = release_pat.match(link)
        if m:
            return link, int(m.group(1)), int(m.group(2))

    # Fallback: scan raw HTML text (handles non-linked entries)
    m = re.search(r"RELEASE_DFAM_(\d+)_(\d+)", html)
    if m:
        return m.group(0), int(m.group(1)), int(m.group(2))

    raise RuntimeError(f"No RELEASE_DFAM tag found at {release_root_url}")


def fetch_file_listing(base_url):
    """Return list of .h5.gz and .h5.gz.md5 filenames from the directory listing."""
    html = _fetch_html(base_url)
    parser = _LinkParser()
    parser.feed(html)
    return [f for f in parser.links if ".h5.gz" in f]


def parse_components(files):
    """
    Parse filenames into {component: [partition, ...]} and detect the release prefix.

    File patterns:
        {prefix}.{partition}.h5.gz              -- root
        {prefix}.{component}.{partition}.h5.gz  -- named component
    """
    gz_files = [f for f in files if f.endswith(".h5.gz")]
    root_pat = re.compile(r"^(dfam\w+?)\.(\d+)\.h5\.gz$")
    comp_pat = re.compile(r"^(dfam\w+?)\.((?:(?:un)?curated)\.(?:consensus|hmm))\.(\d+)\.h5\.gz$")

    prefix = None
    components = defaultdict(set)

    for f in gz_files:
        m = comp_pat.match(f)
        if m:
            prefix = m.group(1)
            components[m.group(2)].add(int(m.group(3)))
            continue
        m = root_pat.match(f)
        if m:
            prefix = m.group(1)
            components["root"].add(int(m.group(2)))

    return prefix, {k: sorted(v) for k, v in components.items()}


def _parse_range_selection(text, available_set, comp):
    """Parse 'all', '0', '0-5', '0,2,5-10' against a set of available ints."""
    text = text.strip().lower()
    if text == "all":
        return sorted(available_set)

    result = set()
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            parts = token.split("-", 1)
            try:
                lower_bound, higher_bound = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            for i in range(lower_bound, higher_bound + 1):
                if i not in available_set:
                    raise ValueError(f"Partition {i} not available for '{comp}'. Should be one of: {available_set}")
            result.update(p for p in available_set if lower_bound <= p <= higher_bound)
        else:
            try:
                p = int(token)
            except ValueError:
                raise ValueError(f"Invalid partition value: {token}. Partition values should be integers.")
            if p in available_set:
                result.add(p)
            else:
                raise ValueError(f"Partition {p} not available for '{comp}'. Should be one of: {available_set}")
    return sorted(result)


def get_components(components, component_selection: str, sizes=None, release=None):
    """Component selection; returns list of chosen component names."""
    ordered = [c for c in COMPONENT_ORDER if c in components]
    ordered += [c for c in sorted(components) if c not in ordered]

    print("\nAvailable components:")
    for i, comp in enumerate(ordered, 1):
        parts = components[comp]
        label = COMPONENT_LABELS.get(comp, comp)
        n = len(parts)
        part_range = f"partition {parts[0]}" if n == 1 else f"partitions 0-{parts[-1]} ({n} total)"
        hint = _size_hint(sizes.get(comp) if sizes else None, n)
        req = " [required]" if comp == "root" else " [optional]"
        print(f"  {i:2}. {label}{req}")
        print(f"        [{comp}] -- {part_range}{hint}")

    if sizes:
        total = 0
        for comp, parts in components.items():
            entry = sizes.get(comp)
            if entry:
                first, last = entry
                avg = first if last is None else (first + last) / 2
                total += avg * len(parts)
        if total:
            label = f"Dfam {release}" if release else "complete"
            print(f"\n  Complete {label} download would be ~{_fmt_size(total)} compressed.")

    #print("\nEnter numbers to download (e.g. '1,3' or 'all'):")
    #print("  (Components with multiple partitions will prompt for partition selection.)")
    print(f"Selected components: {component_selection}")
    if component_selection == "all":
        return ordered
    try:
        chosen = []
        for tok in component_selection.split(","):
            idx = int(tok.strip()) - 1
            if 0 <= idx < len(ordered):
                chosen.append(ordered[idx])
            else:
                raise IndexError(f"Invalid partition selection: {tok}. Partition values must be between 1 and {len(ordered)}")
        if not chosen:
            raise ValueError("No partition selected.")
        return chosen
    except (ValueError, IndexError):
        raise ValueError("Invalid selection")


def _prompt_partitions(comp, available, selected_partitions: str, size_entry=None):
    """Interactive partition selection for a single component."""
    label = COMPONENT_LABELS.get(comp, comp)
    n = len(available)
    hint = _size_hint(size_entry, n)

    if n == 1:
        print(f"\n  {label}: only 1 partition ({available[0]}), selecting it.")
        return available

    #print(f"\n  {label} has {n} partitions (0-{available[-1]}).{hint}")
    #print("  Enter partitions to download: 'all', '0', '0-5', '0,2,5-10', etc.")
    avail_set = set(available)
    result = _parse_range_selection(selected_partitions, avail_set, comp)
    if not result:
        raise ValueError("Invalid selection")
    return result


def _md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _run(cmd: list[str], shell: bool = False):
    subprocess.run(cmd, shell=shell, check=True)


def _download(url: str, data_path: str, ncpus: int):
    cmd = [
        'aria2c', 
        '-s', str(ncpus), 
        '-x', str(ncpus),  
        '--optimize-concurrent-downloads', 
        '--check-integrity=true', 
        '--dir', 
        str(data_path), 
        url
    ]
    _run(cmd)


def _decompress(file: str):
    cmd = ['pigz', '-df', str(file)]
    _run(cmd)


def _head_size(url):
    """Return Content-Length for url via HEAD, or None on failure."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=10) as resp:
            val = resp.headers.get("Content-Length")
            return int(val) if val else None
    except Exception:
        return None


def _fmt_size(n):
    """Format bytes as a human-readable string."""
    if n >= 1_073_741_824:
        return f"{n / 1_073_741_824:.1f} GB"
    if n >= 1_048_576:
        return f"{n / 1_048_576:.0f} MB"
    return f"{n / 1024:.0f} KB"


def fetch_component_sizes(base_url, prefix, components):
    """
    Return {comp: (first_size, last_size)} using HEAD requests.
    For single-partition components last_size is None.
    Missing entries mean all requests failed.
    """
    print("Fetching partition size estimates...", flush=True)
    sizes = {}
    for comp, parts in components.items():
        def _gz(p, _comp=comp):
            return (f"{prefix}.{p}.h5.gz" if _comp == "root"
                    else f"{prefix}.{_comp}.{p}.h5.gz")
        first = _head_size(base_url + _gz(parts[0]))
        if first is None:
            continue
        last = _head_size(base_url + _gz(parts[-1])) if len(parts) > 1 else None
        sizes[comp] = (first, last)
    return sizes


def _size_hint(entry, n_parts):
    """Format a size hint string from a (first_sz, last_sz) entry."""
    if not entry:
        return ""
    first, last = entry
    if n_parts == 1:
        return f"  ({_fmt_size(first)})"
    if last is None:
        return f"  (~{_fmt_size(first)} compressed per partition)"
    lo, hi = min(first, last), max(first, last)
    if lo == hi:
        return f"  (~{_fmt_size(lo)} compressed per partition)"
    return f"  ({_fmt_size(lo)} to {_fmt_size(hi)} compressed per partition)"


def _read_md5_file(path):
    """Read the first whitespace-delimited token from an md5 file."""
    with open(path) as f:
        return f.read().split()[0].lower()


def process_one(base_url, out_dir, prefix, component, partition, ncpus: int, dry_run=False):
    """
    Download, validate, and decompress one partition file.
    Returns True on success (or already done), False on error.
    """
    stem = (
        f"{prefix}.{partition}.h5"
        if component == "root"
        else f"{prefix}.{component}.{partition}.h5"
    )
    gz_name  = stem + ".gz"
    md5_name = gz_name + ".md5"

    final = os.path.join(out_dir, stem)
    gz    = os.path.join(out_dir, gz_name)
    md5f  = os.path.join(out_dir, md5_name)

    tag = f"{component}/partition-{partition}"

    # Already fully decompressed?
    if os.path.exists(final):
        print(f"  {tag}: already present ({stem}), skipping.")
        return True

    if dry_run:
        print(f"  [dry-run] {gz_name}")
        return True

    # Fetch md5 sidecar
    if not os.path.exists(md5f):
        try:
            _download(base_url + md5_name, out_dir, ncpus)
        except Exception as e:
            print(f"  ERROR fetching {md5_name}: {e}")
            return False

    try:
        expected = _read_md5_file(md5f)
    except Exception as e:
        print(f"  ERROR reading {md5f}: {e}")
        return False

    # Download gz if missing or corrupt
    need_dl = True
    if os.path.exists(gz):
        print(f"    Verifying existing {gz_name}...", end=" ", flush=True)
        if _md5sum(gz) == expected:
            print("OK")
            need_dl = False
        else:
            print("MISMATCH - re-downloading")
            os.remove(gz)

    if need_dl:
        try:
            _download(base_url + gz_name, out_dir, ncpus)
        except Exception as e:
            print(f"  ERROR downloading {gz_name}: {e}")
            return False

        actual = _md5sum(gz)
        if actual != expected:
            print(f"  ERROR: MD5 mismatch for {gz_name}")
            print(f"    expected: {expected}")
            print(f"    actual:   {actual}")
            os.remove(gz)
            return False
        print(f"    MD5 OK: {gz_name}")
        
    # Decompress
    print(f"    Decompressing {gz_name} ...", end=" ", flush=True)
    _decompress(gz)

    return True



def parse_famdb_output(file: str):
    components_to_partitions = {}
    with open(file) as fin:
        component = None
        while True:
            try:
                line = next(fin).strip('\n')
            except StopIteration:
                break
            # lines of interest start with a space
            if not line.startswith(" "):
                continue
            # getting words in the line
            words = [s for s in line.split(' ') if s != '']
            # if first word is not partition, it is the first line of a component
            if words[0] != 'partition':
                # getting component
                component_words = []
                for i, word in enumerate(words):
                    if word == 'partition':
                        break
                    component_words.append(word)
                component = FAMDB_OUTPUT_KEYS_TO_COMPONENTS[' '.join(component_words)]
                # the value of the partition is the word right after 'partition'
                components_to_partitions[component] = [words[i + 1]]
            else:
                if component is None:
                    raise ValueError(f"Unknown component as line {line}")
                # the value of the partition is the word right after 'partition'
                components_to_partitions[component].append(words[1])

    # turning lists into strings
    return  {k: ','.join(v) for k, v in components_to_partitions.items()}
    


def main():
    ap = argparse.ArgumentParser(
        description="Download, verify, and decompress FamDB component files from Dfam."
    )
    ap.add_argument(
        "-o", "--output-dir",
        default=None,
        required=True,
        metavar="DIR",
        help="Destination directory (default: Libraries/famdb or famdb.conf setting)",
    )
    ap_root = ap.add_mutually_exclusive_group()
    ap_root.add_argument(
        "--fetch-root",
        action="store_true",
        help="Download Dfam root"
    )
    ap_famdb_out = ap.add_mutually_exclusive_group()
    ap_famdb_out.add_argument(
        "--famdb-check-output",
        dest="famdb_check_output_file",
        help="File containing output of famdb check, used for selection of components and partitions"
    )
    ap.add_argument(
        "-u", "--url",
        default=DEFAULT_URL,
        metavar="URL",
        help=f"Base URL of FamDB release directory (default: {DEFAULT_URL})",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without downloading anything",
    )
    ap.add_argument("--ncpus", required=True, type=int, help='Number of CPUs to use for downloading.')
    args = ap.parse_args()

    if not args.fetch_root and not args.famdb_check_output_file:
        raise RuntimeError("Either --fetch-root or --famdb-check-output must be specified")

    base_url = args.url.rstrip("/") + "/"
    print("#")
    print("# download_dfam.py")
    print("#")
    print("#")
    print("# To identify the minimal set of partitions to download for a given")
    print("# species or taxon, simply download the root partition first, and then")
    print("# query the release details using:")
    print("#    ./famdb.py check <species or taxon>")
    print("#")

    # CHANGED: take args.output_dir as it is
    out_dir = args.output_dir

    print(f"Output directory: {out_dir}")

    # Navigate 2 directory levels up from the FamDB URL to find the release root
    release_root = base_url.rstrip("/").rsplit("/", 2)[0] + "/"
    try:
        tag, major, minor = fetch_release_tag(release_root)
    except RuntimeError as e:
        print(f"ERROR: Could not determine release version: {e}")
        sys.exit(1)

    print(f"Release: {tag}")

    MIN_MAJOR, MIN_MINOR = 4, 0
    if (major, minor) < (MIN_MAJOR, MIN_MINOR):
        print(f"ERROR: Release {tag} predates the minimum supported release "
              f"RELEASE_DFAM_{MIN_MAJOR}_{MIN_MINOR}.")
        print("Please report this problem to help@dfam.org")
        sys.exit(1)

    print(f"Fetching file listing from {base_url} ...")
    try:
        files = fetch_file_listing(base_url)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    if not files:
        print("ERROR: No .h5.gz files found at the given URL.")
        sys.exit(1)

    prefix, components = parse_components(files)
    if not prefix:
        print("ERROR: Could not determine release prefix from filenames.")
        sys.exit(1)

    existing_prefix = _detect_existing_prefix(out_dir)
    if existing_prefix and existing_prefix != prefix:
        print(f"\nWARNING: '{out_dir}' already contains files from a different release")
        print(f"         (found: {existing_prefix!r}, downloading: {prefix!r}).")
        print(f"         Mixing releases in the same directory will corrupt the database.")
        raw = input("Continue anyway? [y/N] ").strip().lower()
        if raw not in ("y", "yes"):
            print("Aborted.")
            sys.exit(0)

    sizes = fetch_component_sizes(base_url, prefix, components)

    if args.fetch_root:
        component_selection = "1"
        components_to_partitions = {"root": "0"}
    elif args.famdb_check_output_file:
        components_to_partitions = parse_famdb_output(args.famdb_check_output_file)
        # indices are 0-based in python while we expect 1-bases indices for components
        component_selection = ','.join([str(COMPONENT_ORDER.index(comp) + 1) for comp in components_to_partitions])

    # CHANGED: replaced _prompt_components by get_components
    selected_components = get_components(
        components, 
        component_selection, 
        sizes, 
        release=f"{major}.{minor}"
    )

    work = []
    for comp in selected_components:
        chosen = _prompt_partitions(
            comp, 
            components[comp], 
            components_to_partitions[comp],
            sizes.get(comp) if sizes else None
        )
        for part in chosen:
            work.append((comp, part))

    if not work:
        print("Nothing selected.")
        sys.exit(0)

    root_in_work = any(comp == "root" for comp, _ in work)
    root_on_disk = os.path.exists(os.path.join(out_dir, f"{prefix}.0.h5"))
    if not root_in_work and not root_on_disk:
        print("\nNote: Root partition is required and has been added to the download.")
        work.insert(0, ("root", 0))

    print(f"\n{'[dry-run] ' if args.dry_run else ''}Downloading {len(work)} file(s) to: {out_dir}\n")
    if not args.dry_run:
        os.makedirs(out_dir, exist_ok=True)

    errors = []
    for i, (comp, part) in enumerate(work, 1):
        print(f"[{i}/{len(work)}] {comp} partition {part}")
        ok = process_one(base_url, out_dir, prefix, comp, part, args.ncpus, dry_run=args.dry_run)
        if not ok:
            errors.append((comp, part))

    print("\n--- Summary ---")
    succeeded = len(work) - len(errors)
    print(f"  {succeeded}/{len(work)} completed successfully")
    if errors:
        print("  Failed:")
        for comp, part in errors:
            print(f"    {comp} partition {part}")
        sys.exit(1)


if __name__ == "__main__":
    main()