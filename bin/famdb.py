#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Usage: famdb.py [-h] [-l LOG_LEVEL] [-i DB_DIR] command ...

    Queries or modifies the contents of a famdb file. For more detailed help
    and information about program options, run `famdb.py --help` or
    `famdb.py <command> --help`.

    This program can also be used as a module. It provides classes and methods
    for working with FamDB files, which contain Transposable Element (TE)
    families and associated taxonomy data.

    # Classes
        Family: Metadata and model of a TE family.
        FamDB: HDF5-based format for storing Family objects.

SEE ALSO:
    Dfam: http://www.dfam.org

AUTHOR(S):
    Anthony Gray <anthony.gray@systemsbiology.org>
    Jeb Rosen <jeb.rosen@systemsbiology.org>

LICENSE:
    This code may be used in accordance with the Creative Commons
    Zero ("CC0") public domain dedication:
    https://creativecommons.org/publicdomain/zero/1.0/

DISCLAIMER:
    This software is provided ``AS IS'' and any express or implied
    warranties, including, but not limited to, the implied warranties of
    merchantability and fitness for a particular purpose, are disclaimed.
    In no event shall the authors or the Dfam consortium members be
    liable for any direct, indirect, incidental, special, exemplary, or
    consequential damages (including, but not limited to, procurement of
    substitute goods or services; loss of use, data, or profits; or
    business interruption) however caused and on any theory of liability,
    whether in contract, strict liability, or tort (including negligence
    or otherwise) arising in any way out of the use of this software, even
    if advised of the possibility of such damage.
"""

import argparse
import json
import logging
import os
import re
import sys
import traceback

LOGGER = logging.getLogger(__name__)

from famdb_globals import (
    FILE_DESCRIPTION,
    FAMILY_FORMATS_EPILOG,
    SINGLE_FAMILY_FORMATS_EPILOG,
    MISSING_FILE,
    HELP_URL,
    COMPONENT_CC,
    COMPONENT_CH,
    COMPONENT_UC,
    COMPONENT_UH,
    COMPONENT_TYPES,
    resolve_db_dir,
)
from famdb_classes import FamDB
from famdb_helper_methods import filter_curated


def _format_partition(partition_dict):
    """Format a partition dict as a compact string, e.g. 'cc:0,ch:1'. Returns 'N/A' if None."""
    if partition_dict is None:
        return "N/A"
    parts = [f"{k}:{v}" for k, v in partition_dict.items() if v is not None]
    return ",".join(parts) if parts else "N/A"


# Command-line utilities
def command_info(args):
    """The 'info' command displays some of the stored metadata."""
    args.db_dir.print_info(history=args.history)


def command_names(args):
    """The 'names' command displays all names of all taxa that match the search term."""

    entries = []
    entries += args.db_dir.resolve_names(args.term)

    if args.format == "pretty":
        print(
            "Partition key - index of the component partition containing this taxon's families:\n"
            "  cc = Curated Consensus    ch = Curated HMMs\n"
            "  uc = Uncurated Consensus  uh = Uncurated HMMs\n"
        )
        prev_exact = None
        for tax_id, is_exact, partition, names in entries:
            if is_exact != prev_exact:
                if is_exact:
                    print("Exact Matches\n=============")
                else:
                    if prev_exact:
                        print()
                    print("Non-exact Matches\n=================")
                prev_exact = is_exact

            print(
                f"Taxon: {tax_id}, Partition: {_format_partition(partition)}, Names: {', '.join([f'{n[1]} ({n[0]})' for n in names])}"
            )

    elif args.format == "json":
        obj = []
        for tax_id, is_exact, partition, names in entries:
            names_obj = [{"kind": name[0], "value": name[1]} for name in names]
            obj += [{"id": tax_id, "partition": partition, "names": names_obj}]
        print(json.dumps(obj))
    else:
        raise ValueError("Unimplemented names format: %s" % args.format)


def _taxon_installed_count(file, tax_id, model, curated_only=False, uncurated_only=False):
    """Returns (installed_count, total_count) for families directly assigned to tax_id.
    model must be 'consensus' or 'hmm'."""
    accessions = file.get_families_for_taxon(tax_id, curated_only, uncurated_only)
    total = len(accessions)
    if total == 0:
        return 0, 0

    if model == "hmm":
        cur_ct, uncur_ct = COMPONENT_CH, COMPONENT_UH
    else:
        cur_ct, uncur_ct = COMPONENT_CC, COMPONENT_UC

    cur_part = file.files[0].get_partition_for_taxon(tax_id, cur_ct)
    uncur_part = file.files[0].get_partition_for_taxon(tax_id, uncur_ct)
    cur_installed = cur_part is not None and cur_part in file.components.get(cur_ct, {})
    uncur_installed = uncur_part is not None and uncur_part in file.components.get(uncur_ct, {})

    if curated_only:
        return (total if cur_installed else 0), total
    if uncurated_only:
        return (total if uncur_installed else 0), total

    if cur_installed and uncur_installed:
        return total, total
    if not cur_installed and not uncur_installed:
        return 0, total

    installed = sum(
        1 for acc in accessions
        if (filter_curated(acc, True) and cur_installed)
        or (not filter_curated(acc, True) and uncur_installed)
    )
    return installed, total


def print_lineage_tree(
    file,
    tree,
    gutter_self,
    gutter_children,
    curated_only=False,
    uncurated_only=False,
    model=None,
):
    """Pretty-prints a lineage tree with box drawing characters."""

    if not tree:
        return
    if type(tree) == str:
        tax_id = tree
        children = []
    else:
        tax_id = tree[0]
        children = tree[1:]
    name, _tax_partition = file.get_taxon_name(tax_id, "scientific name")
    if name != "Not Found":
        if model is not None:
            installed, total = _taxon_installed_count(file, tax_id, model, curated_only, uncurated_only)
            count = f"[{installed}/{total}]"
        else:
            total = len(file.get_families_for_taxon(tax_id, curated_only, uncurated_only))
            count = f"[{total}]"
        print(f"{gutter_self}{tax_id} {name} {count}")

    # All but the last child need a downward-pointing line that will link up
    # to the next child, so this is split into two cases
    if len(children) > 1:
        for child in children[:-1]:
            print_lineage_tree(
                file,
                child,
                gutter_children + "├─",
                gutter_children + "│ ",
                curated_only,
                uncurated_only,
                model,
            )

    if children:
        print_lineage_tree(
            file,
            children[-1],
            gutter_children + "└─",
            gutter_children + "  ",
            curated_only,
            uncurated_only,
            model,
        )


def print_lineage_semicolons(
    file,
    tree,
    parent_name,
    starting_at,
    curated_only=False,
    uncurated_only=False,
):
    """
    Prints a lineage tree as a flat list of semicolon-delimited names.

    In order to print the correct lineage string, the available tree must
    be "complete" even if ancestors were not specified to build up the
    string starting from "root". 'starting_at' specifies the first taxa
    (in the descending direction) to actually be output.
    """
    if not tree:
        return

    tax_id = tree[0]
    children = tree[1:]
    name, tax_partition = file.get_taxon_name(tax_id, "scientific name")

    if name != "Not Found":
        if parent_name:
            name = parent_name + ";" + name

        if starting_at == tax_id:
            starting_at = None

        if not starting_at:
            fams = file.get_families_for_taxon(
                tax_id, curated_only, uncurated_only
            )
            count = f"[{len(fams)}]" if fams is not None else "[?]"
            print(f"{tax_id}({_format_partition(tax_partition)}): {name} {count}")

        for child in children:
            print_lineage_semicolons(
                file,
                child,
                name,
                starting_at,
                curated_only,
                uncurated_only,
            )


def get_lineage_totals(
    file,
    tree,
    target_id,
    curated_only=False,
    uncurated_only=False,
    model="consensus",
    seen=None,
    seen_present=None,
):
    """
    Recursively calculates the total number of families
    on ancestors and descendants of 'target_id' in the given 'tree'.

    'seen' is required to track families that are present on multiple
    lineages due to horizontal transfer and ensure each family
    is only counted one time, either as an ancestor or a descendant.

    Returns [ancestor_count, descendant_count, ancestor_present, descendant_present]
    where *_present reflects only families available in locally installed partitions
    for the given model type ('consensus' or 'hmm').
    """
    if not seen:
        seen = set()
    if seen_present is None:
        seen_present = set()

    tax_id = tree[0]
    children = tree[1:]
    accessions = file.get_families_for_taxon(tax_id, curated_only, uncurated_only)

    count_here = 0
    present_here = 0
    if accessions:
        for acc in accessions:
            if acc not in seen:
                seen.add(acc)
                count_here += 1
            if acc not in seen_present:
                is_curated = filter_curated(acc, True)
                if model == "hmm":
                    ct = COMPONENT_CH if is_curated else COMPONENT_UH
                else:
                    ct = COMPONENT_CC if is_curated else COMPONENT_UC
                part_num = file.files[0].get_partition_for_taxon(tax_id, ct)
                if part_num is not None and part_num in file.components.get(ct, {}):
                    seen_present.add(acc)
                    present_here += 1

    if target_id == tax_id:
        target_id = None

    counts = [0, 0, 0, 0]
    for child in children:
        new_counts = get_lineage_totals(
            file,
            child,
            target_id,
            curated_only,
            uncurated_only,
            model,
            seen,
            seen_present,
        )
        counts[0] += new_counts[0]
        counts[1] += new_counts[1]
        counts[2] += new_counts[2]
        counts[3] += new_counts[3]

    if target_id is None:
        counts[1] += count_here
        counts[3] += present_here
    else:
        counts[0] += count_here
        counts[2] += present_here
    return counts


def command_lineage(args):
    """The 'lineage' command outputs ancestors and/or descendants of the given taxon."""

    target_id, partition = args.db_dir.resolve_one_species(args.term)

    if not target_id:
        print(f"No species found for search term '{args.term}'", file=sys.stderr)
        return
    if target_id == "Ambiguous":
        return
    tree = args.db_dir.get_lineage(
        target_id,
        descendants=args.descendants,
        ancestors=args.ancestors or args.format == "semicolon",
        complete=args.complete or args.format == "semicolon",
    )
    if not tree:
        return
    if args.format == "pretty":
        if args.model is not None:
            if args.curated:
                count_note = "curated family consensus sequences"
            elif args.uncurated:
                count_note = "uncurated (DR) family consensus sequences"
            else:
                count_note = "curated (DF) and uncurated (DR) families"
        else:
            if args.curated:
                count_note = "curated (DF) families"
            elif args.uncurated:
                count_note = "uncurated (DR) families"
            else:
                count_note = "curated (DF) and uncurated (DR) families"
        if args.model is not None:
            print(
                f"# Format: <NCBI tax ID> <scientific name> [<# families_installed>/<# families in Dfam {args.db_dir.db_version}>]\n"
                f"#        where counts represent {count_note}\n"
            )
        else:
            print(
                f"# Format: <NCBI tax ID> <scientific name> [<# families>]\n"
                f"#        where counts represent {count_note}\n"
            )
        print_lineage_tree(
            args.db_dir,
            tree,
            "",
            "",
            args.curated,
            args.uncurated,
            args.model,
        )
    elif args.format == "semicolon":
        print_lineage_semicolons(
            args.db_dir, tree, "", target_id, args.curated, args.uncurated
        )
    elif args.format == "totals":
        totals = get_lineage_totals(
            args.db_dir, tree, target_id, args.curated, args.uncurated,
            args.model or "consensus",
        )
        print(
            f"{totals[0]} entries in ancestors; {totals[1]} lineage-specific entries; "
            f"{totals[2]} ancestral entries present; {totals[3]} lineage-specific entries present"
        )
    else:
        raise ValueError("Unimplemented lineage format: %s" % args.format)


def command_check(args):
    """The 'check' command reports which component partitions are locally installed for a taxon."""

    target_id, _ = args.db_dir.resolve_one_species(args.term)
    if not target_id:
        print(f"No species found for search term '{args.term}'", file=sys.stderr)
        return
    if target_id == "Ambiguous":
        return

    tax_name, _ = args.db_dir.get_taxon_name(target_id, "scientific name")
    print(f"\nPartition check for '{tax_name}' (tax id: {target_id}):\n")

    component_labels = {
        COMPONENT_CC: "Curated Consensus",
        COMPONENT_CH: "Curated HMMs",
        COMPONENT_UC: "Uncurated Consensus",
        COMPONENT_UH: "Uncurated HMMs",
    }
    ct_key = {
        COMPONENT_CC: "cc",
        COMPONENT_CH: "ch",
        COMPONENT_UC: "uc",
        COMPONENT_UH: "uh",
    }

    components_to_check = args.component if args.component else COMPONENT_TYPES

    # Collect all partitions needed: one per unique partition number across the
    # full ancestor lineage (ancestors supply families applicable to this taxon too).
    lineage = args.db_dir.get_lineage_path(target_id)
    needed = {ct: set() for ct in components_to_check}
    for _name, part_dict in lineage:
        if part_dict is None:
            continue
        for ct in components_to_check:
            pn = part_dict.get(ct)
            if pn is not None:
                needed[ct].add(pn)

    max_label = max(len(component_labels[ct]) for ct in components_to_check)

    for ct in components_to_check:
        label = component_labels[ct]
        partitions = sorted(needed[ct])
        if not partitions:
            print(f"  {label:<{max_label}}  N/A  (no families for this taxon or its ancestors)")
        else:
            for i, part_num in enumerate(partitions):
                display_label = label if i == 0 else ""
                fm_key = f"{ct_key[ct]}.{part_num}"
                fm_entry = args.db_dir.file_map.get(fm_key, {})
                root_name = fm_entry.get("T_root_name", "")
                part_label = f"partition {part_num} [{root_name}]:" if root_name else f"partition {part_num}:"
                leaf = args.db_dir.components[ct].get(part_num)
                if leaf is not None:
                    print(f"  {display_label:<{max_label}}  {part_label}  present")
                else:
                    filename = fm_entry.get("filename", f"partition {part_num}")
                    print(f"  {display_label:<{max_label}}  {part_label}  MISSING  [{filename}]")
    print()


def print_families(args, families, header, species=None):
    """
    Prints each family in 'families', optionally with a copyright header. The
    format is determined by 'args.format' and additional data (such as
    taxonomy) is taken from 'args.db_dir'.

    If 'species' is provided and the format is "hmm_species", it is the id of
    the taxa whose species-specific thresholds should be substituted into the
    GA, NC, and TC lines of the HMM.
    """

    # These args are only available with the "families" command. When
    # print_families is called by the "family" command, accessing e.g.
    # args.stage directly raises an AttributeError
    # TODO: consider reworking argument passing to avoid this workaround
    add_reverse_complement = getattr(args, "add_reverse_complement", False)
    include_class_in_name = getattr(args, "include_class_in_name", False)
    require_general_threshold = getattr(args, "require_general_threshold", False)
    stage = getattr(args, "stage", None)

    if header:
        db_info = args.db_dir.get_metadata()
        if db_info:
            copyright_text = db_info["copyright"]
            # Add appropriate comment character to the copyright header lines
            if "hmm" in args.format:
                copyright_text = re.sub("(?m)^", "#   ", copyright_text)
            elif "fasta" in args.format:
                copyright_text = None
            elif "embl" in args.format:
                copyright_text = re.sub("(?m)^", "CC   ", copyright_text)
            if copyright_text:
                print(copyright_text)

    for family in families:
        if args.format == "summary":
            if include_class_in_name:
                name = family.name or family.accession
                rm_class = family.repeat_type
                if family.repeat_subtype:
                    rm_class += "/" + family.repeat_subtype
                family.name = name + "#" + rm_class
            entry = str(family) + "\n"
        elif args.format == "hmm":
            entry = family.to_dfam_hmm(
                args.db_dir,
                include_class_in_name=include_class_in_name,
                require_general_threshold=require_general_threshold,
            )
        elif args.format == "hmm_species":
            entry = family.to_dfam_hmm(
                args.db_dir,
                species,
                include_class_in_name=include_class_in_name,
                require_general_threshold=require_general_threshold,
            )
        elif (
            args.format == "fasta"
            or args.format == "fasta_name"
            or args.format == "fasta_acc"
        ):
            use_accession = args.format == "fasta_acc"

            buffers = []
            if stage and family.buffer_stages:
                for spec in family.buffer_stages.split(","):
                    if "[" in spec:
                        matches = re.match(r"(\d+)\[(\d+)-(\d+)\]", spec.strip())
                        if matches:
                            if stage == int(matches.group(1)):
                                buffers += [
                                    [int(matches.group(2)), int(matches.group(3))]
                                ]
                        else:
                            LOGGER.warning(
                                "Ingored invalid buffer specification: '%s'",
                                spec.strip(),
                            )
                    else:
                        buffers += [stage == int(spec)]

            if not buffers:
                buffers += [None]

            entry = ""
            for buffer_spec in buffers:
                entry += (
                    family.to_fasta(
                        args.db_dir,
                        use_accession=use_accession,
                        include_class_in_name=include_class_in_name,
                        buffer=buffer_spec,
                    )
                    or ""
                )

                if add_reverse_complement:
                    entry += (
                        family.to_fasta(
                            args.db_dir,
                            use_accession=use_accession,
                            include_class_in_name=include_class_in_name,
                            do_reverse_complement=True,
                            buffer=buffer_spec,
                        )
                        or ""
                    )
        elif args.format == "embl":
            entry = family.to_embl(args.db_dir)
        elif args.format == "embl_meta":
            entry = family.to_embl(args.db_dir, include_meta=True, include_seq=False)
        elif args.format == "embl_seq":
            entry = family.to_embl(args.db_dir, include_meta=False, include_seq=True)
        else:
            raise ValueError("Unimplemented family format: %s" % args.format)

        if entry:
            print(entry, end="")


def _diagnose_missing_accessions(db_dir, missing_accessions, term):
    """
    Classify accessions that were in the index but not found in any component
    file, then emit a targeted error or warning.

    If the missing accessions are all of one curatedness type and the
    corresponding component files are simply absent, that is an expected
    partial-installation situation and deserves a clear ERROR.  Anything
    else is an unexpected data-integrity problem and gets a WARNING.
    """
    if not missing_accessions:
        return

    uncurated_missing = [a for a in missing_accessions if not filter_curated(a, True)]
    curated_missing   = [a for a in missing_accessions if filter_curated(a, True)]

    uc_present = bool(db_dir.components[COMPONENT_UC]) or bool(db_dir.components[COMPONENT_UH])
    cc_present = bool(db_dir.components[COMPONENT_CC]) or bool(db_dir.components[COMPONENT_CH])

    unexplained = []

    if uncurated_missing and not uc_present:
        ex = uncurated_missing[0]
        LOGGER.error(
            "%d uncurated accession(s) (e.g. %s) were requested for '%s' however "
            "the uncurated component is not present on this system.  Please use "
            "'./famdb.py check \"%s\"' to locate missing partitions or change your "
            "--format/--curated/--uncurated options.",
            len(uncurated_missing), ex, term, term,
        )
    else:
        unexplained.extend(uncurated_missing)

    if curated_missing and not cc_present:
        ex = curated_missing[0]
        LOGGER.error(
            "%d curated accession(s) (e.g. %s) were requested for '%s' however "
            "the curated component is not present on this system.  Please use "
            "'./famdb.py check \"%s\"' to locate missing partitions or change your "
            "--format/--curated/--uncurated options.",
            len(curated_missing), ex, term, term,
        )
    else:
        unexplained.extend(curated_missing)

    for acc in unexplained:
        LOGGER.warning(
            "Accession %s found in index but missing from all loaded component files "
            "(possible data integrity issue)", acc
        )


def command_family(args):
    """The 'family' command outputs a single family by name or accession."""
    family = args.db_dir.get_family_by_accession_merged(args.accession)
    if not family:
        family = args.db_dir.get_family_by_name(args.accession)

    if family:
        print_families(args, [family], False)
    else:
        _diagnose_missing_accessions(args.db_dir, [args.accession], args.accession)


def command_families(args):
    """The 'families' command outputs all families associated with the given taxon."""
    target_id, _ = args.db_dir.resolve_one_species(args.term)
    if not target_id:
        print(f"No species found for search term '{args.term}'", file=sys.stderr)
        return
    elif target_id == "Ambiguous":
        return

    is_hmm = args.format.startswith("hmm")

    # NB: This is speed-inefficient, because get_accessions_filtered needs to
    # read the whole family data even though we read it again right after.
    # However it is *much* more memory-efficient than loading all the family
    # data at once and then sorting by accession.
    accessions = sorted(
        args.db_dir.get_accessions_filtered(
            tax_id=target_id,
            descendants=args.descendants,
            ancestors=args.ancestors,
            curated_only=args.curated,
            uncurated_only=args.uncurated,
            is_hmm=is_hmm,
            stage=args.stage,
            repeat_type=args.repeat_type,
            name=args.name,
        )
    )
    # For formats that need both consensus and pHMM data, use the merged getter
    needs_merge = is_hmm or "embl" in args.format
    getter = args.db_dir.get_family_by_accession_merged if needs_merge else args.db_dir.get_family_by_accession

    missing_accessions = []

    def getter_checked(acc):
        fam = getter(acc)
        if fam is None:
            missing_accessions.append(acc)
        return fam

    families = filter(None, map(getter_checked, accessions))

    header = True if accessions else False
    print_families(args, families, header, target_id)
    _diagnose_missing_accessions(args.db_dir, missing_accessions, args.term)


# RepeatMasker Commands -----------------------------------------------------------------------
def command_fasta_all(args):
    """
    command prints out all curated families in FASTA format
    This command is not documented in the help. It is used to export all of the curated families
    to FASTA format for use by RepeatMasker
    """
    args.format = "fasta_name"
    args.include_class_in_name = True
    print_families(args, args.db_dir.fasta_all("/DF"), True, 1)
    print_families(args, args.db_dir.fasta_all("/Aux"), True, 1)


def command_repeatpeps(args):
    """prints the RepeatPeps file"""
    print(args.db_dir.get_repeatpeps())


def command_edit_description(args):
    """Updates the db description"""
    args.db_dir.update_description(args.new)


def command_append(args):
    """
    The 'append' command reads an EMBL file and appends its entries to an
    existing famdb file.
    """

    lookup = args.db_dir.get_all_taxa_names()
    # infile_lookup = {}
    # with open(args.infile) as file:
    #     infile_lookup = json.load(file)
    # lookup.update(infile_lookup)

    header = None

    def set_header(val):
        nonlocal header
        header = val

    embl_iter = FamDB.read_embl_families(args.infile, lookup, header_cb=set_header)

    message = f"Adding Families From {args.infile.split('/')[-1]}"
    rec = args.db_dir.append_start_changelog(message)

    LOGGER.info(message)
    total_ctr = 0
    added_ctr = 0
    file_counts = {}
    new_val_taxa = set()
    dups = set()
    missing_parts = {}  # {partition_num: count}

    cc_components = args.db_dir.components[COMPONENT_CC]

    for entry in embl_iter:
        # check installation namespace and skip entry if it already exists
        # 2026/02/24: Neglected check against lowercase names
        if entry.accession.lower in args.rb_names or not args.db_dir.check_unique(entry):
            continue

        total_ctr += 1
        acc = entry.accession
        added = False

        # Route each clade to the appropriate CC partition file
        add_leaves = {}   # {partition_num: FamDBLeaf}
        add_taxa = set()
        for clade in entry.clades:
            part_dict = args.db_dir.find_taxon(clade)
            cc_part = part_dict.get("cc") if part_dict else None
            if cc_part is not None and cc_part in cc_components:
                leaf = cc_components[cc_part]
                add_leaves[cc_part] = leaf
                # check if the taxon currently has no families (newly valued)
                if not args.db_dir.get_families_for_taxon(clade):
                    add_taxa.add(clade)
            elif cc_part is not None:
                missing_parts[cc_part] = missing_parts.get(cc_part, 0) + 1

        if not add_leaves:
            LOGGER.debug(f" {acc} not added to local files, no CC partition file found")

        for part_num, leaf in add_leaves.items():
            try:
                leaf.add_family(entry)
                # Update root Lookup/ByTaxon for this family
                args.db_dir.files[0]._add_family_taxon_links(acc, entry.clades)
                LOGGER.debug(f"Added {acc} to CC partition {part_num}")
                if not added:
                    added_ctr += 1
                    added = True
                file_counts[part_num] = file_counts.get(part_num, 0) + 1
            except Exception as e:
                LOGGER.debug(f" Ignoring duplicate entry {entry.accession}: {e}")
                dups.add(entry.accession)

        # track formerly empty clades with new additions
        if added:
            new_val_taxa.update(add_taxa)

    args.db_dir.append_finish_changelog(message, rec)
    args.db_dir.update_changelog(added_ctr, total_ctr, file_counts, args.infile)

    LOGGER.info(f"Added {added_ctr}/{total_ctr} families")
    if dups:
        LOGGER.debug(f" {len(dups)} Duplicate Accesisons: {dups}")
    if missing_parts:
        for part_num in missing_parts:
            LOGGER.info(
                f"FamDB CC Partition {part_num} Not Found. {missing_parts[part_num]} Entries Were Not Included"
            )

    db_info = args.db_dir.get_metadata()

    if args.name:
        db_info["name"] = args.name
    if args.description:
        db_info["description"] += "\n" + args.description

    if header:
        db_info["copyright"] += f"\n\n{header}"

    args.db_dir.set_db_info(
        db_info["name"],
        db_info["db_version"],
        db_info["date"],
        db_info["description"],
        db_info["copyright"],
    )

    # Write the updated counts and metadata
    if new_val_taxa:
        LOGGER.info("Rebuilding Sparse Taxonomy Tree")
        args.db_dir.rebuild_pruned_tree(new_val_taxa)

    LOGGER.info("Finalizing Files")
    args.db_dir.finalize()


def build_args():
    """builds and parses the command line args"""
    parser = argparse.ArgumentParser(
        description=FILE_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-l", "--log_level", default="INFO")

    parser.add_argument("-i", "--db_dir", help="specifies the directory to query")

    parser.add_argument(
        "-e",
        "--exclude_files",
        help="exclude specific files, in the form -e X, or -e X,Y,Z for multiple files. 0 has no effect",
    )

    subparsers = parser.add_subparsers(
        description="""Specifies the kind of query to perform.
For more information on all the possible options for a command, add the --help option after it:
famdb.py families --help
""",
        #  metavar, if specified overrides what shows up on the help line as valid
        #  subcommands.  All subcommands will however be printed in the error message
        #  if a bad subcommand is entered as a possibility, so it doesn't hide it
        #  completely.  This is added to hide the new fasta_all command.
        metavar="{info,names,lineage,check,families,family,append}",
    )
    # INFO --------------------------------------------------------------------------------------------------------------------------------
    p_info = subparsers.add_parser(
        "info", description="List general information about the file."
    )
    p_info.add_argument(
        "--history",
        action="store_true",
        help="List the file changelog in addition to general information",
    )
    p_info.set_defaults(func=command_info)

    # NAMES --------------------------------------------------------------------------------------------------------------------------------
    p_names = subparsers.add_parser(
        "names", description="List the names and taxonomy identifiers of a clade."
    )
    p_names.add_argument(
        "-f",
        "--format",
        default="pretty",
        choices=["pretty", "json"],
        metavar="<format>",
        help="choose output format. The default is 'pretty'. 'json' is more appropriate for scripts.",
    )
    p_names.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or part of a scientific or common name",
    )
    p_names.set_defaults(func=command_names)

    # LINEAGE --------------------------------------------------------------------------------------------------------------------------------
    p_lineage = subparsers.add_parser(
        "lineage",
        description="List the taxonomy tree including counts of families at each clade.",
    )
    p_lineage.add_argument(
        "-a",
        "--ancestors",
        action="store_true",
        help="include all ancestors of the given clade",
    )
    p_lineage.add_argument(
        "-d",
        "--descendants",
        action="store_true",
        help="include all descendants of the given clade",
    )
    p_lineage.add_argument(
        "-k",
        "--complete",
        action="store_true",
        help="include output of taxa without families",
        default=False,
    )
    p_lineage.add_argument(
        "-c",
        "--curated",
        action="store_true",
        help="only tabulate curated families ('DF' records)",
    )
    p_lineage.add_argument(
        "-u",
        "--uncurated",
        action="store_true",
        help="only tabulate uncurated families ('DR' records)",
    )
    p_lineage.add_argument(
        "-f",
        "--format",
        default="pretty",
        choices=["pretty", "semicolon", "totals"],
        metavar="<format>",
        help="choose output format. The default is 'pretty'. 'semicolon' is more appropriate for scripts. 'totals' displays the number of ancestral and lineage-specific families found.",
    )
    p_lineage.add_argument(
        "--model",
        default=None,
        choices=["consensus", "hmm"],
        metavar="<model>",
        help="model type ('consensus' or 'hmm'). In 'pretty' format, enables [present/total] counts showing locally installed families. In 'totals' format, selects which model type to check (defaults to 'consensus').",
    )
    p_lineage.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name",
    )
    p_lineage.set_defaults(func=command_lineage)

    # CHECK --------------------------------------------------------------------------------------------------------------------------------
    p_check = subparsers.add_parser(
        "check",
        description="Check which component partitions are locally installed for a given taxon.",
    )
    p_check.add_argument(
        "--component",
        action="append",
        choices=COMPONENT_TYPES,
        metavar="<component>",
        dest="component",
        help=f"restrict check to one or more component types ({', '.join(COMPONENT_TYPES)}); may be repeated",
    )
    p_check.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name",
    )
    p_check.set_defaults(func=command_check)

    # FAMILIES --------------------------------------------------------------------------------------------------------------------------------
    family_formats = [
        "summary",
        "hmm",
        "hmm_species",
        "fasta_name",
        "fasta_acc",
        "embl",
        "embl_meta",
        "embl_seq",
    ]
    single_family_formats = [f for f in family_formats if f != "hmm_species"]
    family_formats_epilog = FAMILY_FORMATS_EPILOG
    single_family_formats_epilog = SINGLE_FAMILY_FORMATS_EPILOG

    p_families = subparsers.add_parser(
        "families",
        description="Retrieve the families associated \
with a given clade, optionally filtered by additional criteria",
        epilog=family_formats_epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_families.add_argument(
        "-a",
        "--ancestors",
        action="store_true",
        help="include all ancestors of the given clade",
    )
    p_families.add_argument(
        "-d",
        "--descendants",
        action="store_true",
        help="include all descendants of the given clade",
    )
    p_families.add_argument(
        "--stage",
        type=int,
        help="include only families that should be searched in the given stage",
    )
    p_families.add_argument(
        "--class",
        dest="repeat_type",
        type=str,
        help="include only families that have the specified repeat Type/SubType",
    )
    p_families.add_argument(
        "--name",
        type=str,
        help="include only families whose name begins with this search term",
    )
    p_families.add_argument(
        "-u",
        "--uncurated",
        action="store_true",
        help="include only 'uncurated' families (i.e. named DRXXXXXXXXX)",
    )
    p_families.add_argument(
        "-c",
        "--curated",
        action="store_true",
        help="include only 'curated' families (i.e. not named DFXXXXXXXXX)",
    )
    p_families.add_argument(
        "-f",
        "--format",
        default="summary",
        choices=family_formats,
        metavar="<format>",
        help="choose output format.",
    )
    p_families.add_argument(
        "--add-reverse-complement",
        action="store_true",
        help="include a reverse-complemented copy of each matching family; only suppported for fasta formats",
    )
    p_families.add_argument(
        "--include-class-in-name",
        action="store_true",
        help="include the RepeatMasker type/subtype after the name (e.g. HERV16#LTR/ERVL); only supported for hmm and fasta formats",
    )
    p_families.add_argument(
        "--require-general-threshold",
        action="store_true",
        help="skip families missing general thresholds (and log their accessions at the debug log level)",
    )
    p_families.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name",
    )
    p_families.set_defaults(func=command_families)

    # FAMILY --------------------------------------------------------------------------------------------------------------------------------
    p_family = subparsers.add_parser(
        "family",
        description="Retrieve details of a single family.",
        epilog=single_family_formats_epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_family.add_argument(
        "-f",
        "--format",
        default="summary",
        choices=single_family_formats,
        metavar="<format>",
        help="choose output format.",
    )
    p_family.add_argument(
        "accession", help="the accession of the family to be retrieved"
    )
    p_family.set_defaults(func=command_family)

    # APPEND --------------------------------------------------------------------------------------------------------------------------------
    p_append = subparsers.add_parser("append")
    p_append.add_argument("infile", help="the name of the input file to be appended")
    p_append.add_argument(
        "exclusion_file",
        help="the name of the file listing family names to be excluded",
    )
    p_append.add_argument(
        "--name", help="new name for the database (replaces the existing name)"
    )
    p_append.add_argument(
        "--description",
        help="additional database description (added to the existing description)",
    )
    p_append.set_defaults(func=command_append)

    # FASTA ALL --------------------------------------------------------------------------------------------------------------------------------
    p_fasta = subparsers.add_parser("fasta_all")
    p_fasta.set_defaults(func=command_fasta_all)

    # RepeatPeps -------------------------------------------------------------------------------------------------------------------------------
    p_rp = subparsers.add_parser("repeat_peps")
    p_rp.set_defaults(func=command_repeatpeps)

    # Edit Description -------------------------------------------------------------------------------------------------------------------------------
    p_desc = subparsers.add_parser("edit_description")
    p_desc.add_argument("new")
    p_desc.set_defaults(func=command_edit_description)

    return parser


def main():  # ================================================================================================================================
    """Parses command-line arguments and runs the requested command."""

    logging.basicConfig()

    parser = build_args()
    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    write_commands = [command_append, command_edit_description]
    if "func" in args and args.func in write_commands:
        mode = "r+"
    else:
        mode = "r"

    if "term" in args:
        args.term = " ".join(args.term)

    args.db_dir = resolve_db_dir(args.db_dir)

    if not (args.db_dir and os.path.isdir(args.db_dir)):
        LOGGER.error(
            "FamDB data directory not found. Use -i to specify the directory containing "
            "the *.h5 files, or set FAMDB_DATA_DIR in famdb.conf."
        )
        exit(1)

    if hasattr(args, "func") and args.func.__name__ == "command_append":
        if os.path.exists(args.exclusion_file):
            try:
                with open(args.exclusion_file) as f:
                    args.rb_names = set(name.strip() for name in f.readlines())
            except Exception:
                LOGGER.error(f"{args.exclusion_file} could not be parsed.")
                exit(1)
        else:
            LOGGER.error(f"{args.exclusion_file} not found.")
            exit(1)

    try:
        exclude = (
            [n.strip() for n in args.exclude_files.split(",")]
            if args.exclude_files
            else []
        )
        args.db_dir = FamDB(args.db_dir, mode, exclude)
    except:
        args.db_dir = None
        raise

    if not args.db_dir:
        return

    if "func" in args:
        try:
            args.func(args)
        except Exception:
            traceback.print_exc()
    else:
        parser.print_help()


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # This workaround is from
        # https://docs.python.org/3/library/signal.html#note-on-sigpipe

        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE  