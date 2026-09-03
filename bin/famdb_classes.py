import datetime
import gzip
import time
import os
import json
import sys
import re
import h5py
import numpy

from famdb_helper_classes import Family, TaxNode
import logging
LOGGER = logging.getLogger(__name__)

from famdb_globals import (
    FAMDB_VERSION,
    GROUP_FAMILIES,
    GROUP_LOOKUP_BYNAME,
    GROUP_LOOKUP_BYSTAGE,
    GROUP_LOOKUP_BYTAXON,
    GROUP_NODES,
    GROUP_FILE_HISTORY,
    GROUP_REPEATPEPS,
    DATA_CHILDREN,
    DATA_PARENT,
    DATA_VAL_CHILDREN,
    DATA_VAL_PARENT,
    DATA_TAXANAMES,
    DATA_NAMES_CACHE,
    DATA_PARTITION_CACHE,
    COMPONENT_CC,
    COMPONENT_CH,
    COMPONENT_UC,
    COMPONENT_UH,
    COMPONENT_TYPES,
    COMPONENT_META,
    FAMDB_ROOT_FILE_RE,
    FAMDB_COMPONENT_FILE_RE,
    META_DB_VERSION,
    META_DB_DESCRIPTION,
    META_DB_COPYRIGHT,
    META_DB_DATE,
    META_DB_NAME,
    META_CREATED,
    META_META,
    META_UUID,
    META_FILE_INFO,
    META_FAMDB_VERSION,
    META_FILE_MAP,
    DESCRIPTION,
)
from famdb_helper_methods import (
    sanitize_name,
    sounds_like,
    families_iterator,
    filter_curated,
    filter_repeat_type,
    filter_search_stages,
    filter_defined_search_stages,
    filter_name,
    get_family,
    accession_bin,
    is_fasta,
)


class FamDBLeaf:
    """Transposable Element Family and taxonomy database."""

    dtype_str = h5py.special_dtype(vlen=str)

    def __init__(self, filename, mode="r", component_type=None):
        if mode == "r":
            reading = True

            # If we definitely will not be writing to the file, optimistically assume
            # nobody else is writing to it and disable file locking. File locking can
            # be a bit flaky, especially on NFS, and is unnecessary unless there is
            # a parallel writer (which is unlikely for famdb files).
            os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        elif mode == "r+":
            reading = True
        elif mode == "w":
            reading = False
        else:
            raise ValueError(
                f"Invalid file mode. Expected 'r' or 'r+' or 'w', got '{mode}'"
            )

        self.filename = filename
        self.file = h5py.File(filename, mode)
        self.mode = mode

        if (reading and not self.file.attrs.get(META_FAMDB_VERSION)) or (
            reading and not self.version_match()
        ):
            LOGGER.error(
                f"\t Partition {self.get_partition_num()}:This file cannot be read by this version of famdb.py.\n"
                f" Export File Version: {self.file.attrs.get(META_FAMDB_VERSION, 'Not Found')}\n"
                f" FamDB Script Version: {FAMDB_VERSION}\n"
            )
            sys.exit(1)

        if self.mode == "w":
            self.seen = {}
            self.added = {"consensus": 0, "hmm": 0}
            if component_type is not None:
                if component_type not in COMPONENT_TYPES:
                    raise ValueError(f"Invalid component_type '{component_type}'. Expected one of {COMPONENT_TYPES}")
                self.component_type = component_type
                self.file.attrs["component_type"] = component_type
            else:
                self.component_type = None
            self.__write_metadata()
            # ensure lookups exist to avoid random breaking depending on the export data
            self.file.require_group(GROUP_LOOKUP_BYNAME)
            self.file.require_group(GROUP_LOOKUP_BYSTAGE)
        elif self.mode == "r+":
            self.added = self.get_counts()
            self.component_type = self.file.attrs.get("component_type", None)
        else:
            self.component_type = self.file.attrs.get("component_type", None)

    def version_match(self):
        file_version = self.file.attrs.get(META_FAMDB_VERSION)
        file_splits = file_version.split(".")
        file_major = file_splits[0] if file_splits else None

        script_splits = FAMDB_VERSION.split(".")
        script_major = script_splits[0] if script_splits else None

        same_major = file_major == script_major

        if not same_major:
            return False
        return True

    def update_changelog(self, message, verified=False):
        """
        Creates a OtherData/FileHistory/Timestamp/Message/bool
        to record file changes. Defaults to False to show that change is not complete
        """
        time_stamp = str(datetime.datetime.now())
        group = self.file.require_group(GROUP_FILE_HISTORY).require_group(time_stamp)
        group.create_dataset(message, data=numpy.array([verified]))
        return time_stamp

    def _verify_change(self, time_stamp, message):
        """
        Sets the data of a log entry to True, indicating that it was successful
        """
        self.file[GROUP_FILE_HISTORY][time_stamp][message][0] = True

    def _change_logger(func):
        """
        A wrapper method to update and verify the changelog for common methods
        """
        func_to_note = {
            "__write_metadata": "File Initialized",
            "set_metadata": "Metadata Set",
            "add_family": "Family Added",
            "write_repeatpeps": "RepeatPeps Written",
            "write_full_taxonomy": "Taxonomy Nodes Written",
            "write_lookup_bytaxon": "ByTaxon Lookup Written",
            "write_partition_cache": "Partition Cache Written",
            "update_description": "File Description Updated",
        }
        message = func_to_note[func.__name__]

        def wrapper(self, *args, **kwargs):
            time_stamp = self.update_changelog(message)
            func(self, *args, **kwargs)
            self._verify_change(time_stamp, message)

        return wrapper

    # Export Setters ----------------------------------------------------------------------------------------------------
    @_change_logger
    def __write_metadata(self):
        """Sets file data during writing. Called during file creation"""
        self.file.attrs[META_FAMDB_VERSION] = FAMDB_VERSION
        self.file.attrs[META_CREATED] = str(datetime.datetime.now())
        self.file.attrs[META_DB_DESCRIPTION] = DESCRIPTION

    @_change_logger
    def set_metadata(self, partition_num, map_str, name, version, date, copyright_text, is_root=False):
        """
        Sets database metadata for the current file.
        Stores information about other files as json string.
        Sets partition number (key to file info) and bool if is root file or not.
        'partition_num' may be an int (0 for root) or a string like "cc.0", "ch.1".
        """
        self.file.attrs[META_DB_NAME] = name
        self.file.attrs[META_DB_VERSION] = version
        self.file.attrs[META_DB_DATE] = date
        self.file.attrs[META_DB_COPYRIGHT] = copyright_text

        self.file.attrs[META_FILE_INFO] = json.dumps(map_str)

        self.file.attrs["partition_num"] = str(partition_num)
        self.file.attrs["root"] = is_root or partition_num == "0" or partition_num == 0

    def finalize(self):
        """Writes some collected metadata, such as counts, to the database"""
        self.file.attrs["count_consensus"] = self.added["consensus"]
        self.file.attrs["count_hmm"] = self.added["hmm"]

    @_change_logger
    def update_description(self, new_desc):
        """Updates the description. Available to the user and during the append command"""
        self.file.attrs[META_DB_DESCRIPTION] = new_desc

    # Attribute Getters -----------------------------------------------------------------------------------------------
    def get_partition_num(self):
        """Partition num is used as the key in file_info"""
        return self.file.attrs["partition_num"]

    def get_file_info(self):
        """returns dictionary containing information regarding other related files"""
        return json.loads(self.file.attrs[META_FILE_INFO])

    def is_root(self):
        """Tests if file is root file"""
        return self.file.attrs["root"]

    def get_metadata(self):
        """
        Gets file metadata for the current file as a dict with keys
        'famdb_version', 'created', 'partition_name', 'partition_detail',
        'db_name', 'db_version', 'db_date', 'db_description', 'db_copyright'
        """
        if "db_name" not in self.file.attrs:
            return None
        num = self.get_partition_num()
        partition = self.get_file_info()[META_FILE_MAP][str(num)]
        return {
            "famdb_version": self.file.attrs[META_FAMDB_VERSION],
            "created": self.file.attrs[META_CREATED],
            "partition_name": partition["T_root_name"],
            "partition_detail": ", ".join(partition["F_roots_names"]),
            "name": self.file.attrs[META_DB_NAME],
            "db_version": self.file.attrs[META_DB_VERSION],
            "date": self.file.attrs[META_DB_DATE],
            "description": self.file.attrs[META_DB_DESCRIPTION],
            "copyright": self.file.attrs[META_DB_COPYRIGHT],
        }

    def get_history(self):
        """
        Retrieves and concatenates the changelog into a string
        """
        history = self.file.get(GROUP_FILE_HISTORY)
        messages = {stamp: list(history[stamp].keys())[0] for stamp in history.keys()}
        hist_str = f"\n File {self.get_partition_num()}\n"
        for entry in messages:
            hist_str += f"{entry} - {messages[entry]}\n"
        return hist_str

    def get_counts(self):
        """
        Gets counts of entries in the current file as a dict
        with 'consensus', 'hmm'
        """
        return {
            "consensus": self.file.attrs["count_consensus"],
            "hmm": self.file.attrs["count_hmm"],
        }

    # File Utils
    def interrupt_check(self):
        """
        Changelogs Start as False and are flipped to True when complete
        Returns bool if any changes are not confirmed
        """
        interrupted = False
        history = self.file.get(GROUP_FILE_HISTORY)
        for el in history:
            item = history.get(el)
            note = list(item.keys())[0]
            val = item[note][()][0]
            if not val:
                interrupted = True
                break
        return interrupted

    def close(self):
        """Closes this FamDB instance, making further use invalid."""
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    # Data Writing Methods ---------------------------------------------------------------------------------------------
    # Family Methods
    def check_unique(self, family):
        """Verifies that 'family' is uniquely identified by its value of 'key'."""

        # This is awkward. The EMBL files being appended may only have an
        # "accession", but that accession may match the *name* of a family
        # already in Dfam. The accession may also match a family already in
        # Dfam, but with a "v" added.
        # This has been spot-checked and seems to avoid conflicts - Anthony 11/5/24

        # check by accession first
        accession = family.accession
        binned_acc = accession_bin(accession)
        binned_v = accession_bin(accession + "v")

        if self.file.get(f"{binned_acc}/{accession}") or self.file.get(
            f"{binned_v}/{accession}v"
        ):
            return False

        if self.file.get(f"{GROUP_LOOKUP_BYNAME}/{accession}") or self.file.get(
            f"{GROUP_LOOKUP_BYNAME}/{accession}v"
        ):
            return False

        # check for unique name
        if family.name:
            # name_lookup = f"{GROUP_LOOKUP_BYNAME}/{family.name.lower()}" # TODO add case insensitivity
            name_lookup = f"{GROUP_LOOKUP_BYNAME}/{family.name}"
            if self.file.get(name_lookup) or self.file.get(name_lookup + "v"):
                return False

        return True

    # no @_change_logger here to avoid 1000s of history logs. it is called in the methods that call add_family
    def add_family(self, family):
        """Adds the family described by 'family' to the database.

        Fields not relevant to this file's component_type are stripped before
        writing so that consensus files never store pHMM data and vice versa.
        The original Family object is not modified.
        """
        # Verify uniqueness of name and accession.
        # This is important because of the links created to them later.
        if not self.check_unique(family):
            raise Exception(
                f"Family is not unique! Already seen {family.accession} {f'({family.name})' if family.name else ''}"
            )

        # Strip fields that don't belong to this component type.
        # Work on attribute values directly rather than mutating the Family object.
        is_hmm = self.component_type in (COMPONENT_CH, COMPONENT_UH) if self.component_type else False
        is_consensus = self.component_type in (COMPONENT_CC, COMPONENT_UC) if self.component_type else True

        consensus_val = family.consensus if is_consensus or not self.component_type else None
        model_val = family.model if is_hmm or not self.component_type else None

        # Increment counts
        if consensus_val:
            self.added["consensus"] += 1
        if model_val:
            self.added["hmm"] += 1

        # Create the family data
        # In v0.5 we bin the datasets into subgroups to improve performance
        group_path = accession_bin(family.accession)
        dset = self.file.require_group(group_path).create_dataset(
            family.accession, (0,)
        )

        # Set the family attributes, honouring component-type field restrictions.
        # The model is stored separately as a gzip-compressed sibling dataset
        # (not as an attr) so it never needs to be decompressed on load unless
        # the caller explicitly requests it.
        hmm_only_fields = {"model", "max_length", "is_model_masked", "seed_count",
                           "build_method", "search_method", "taxa_thresholds", "general_cutoff"}
        consensus_only_fields = {"consensus"}

        for k in Family.META_LOOKUP:
            if k == "model":
                continue  # model is stored as a compressed sibling dataset below
            if k in hmm_only_fields and is_consensus and self.component_type:
                continue
            if k in consensus_only_fields and is_hmm and self.component_type:
                continue
            value = getattr(family, k)
            if value:
                dset.attrs[k] = value

        # Store HMM model as a gzip-compressed uint8 dataset alongside the
        # family dataset.  family.model may already be raw gzip bytes (from
        # the export pickle path) or a plain string (from EMBL/HMM file import).
        if model_val is not None:
            if isinstance(model_val, bytes):
                compressed = model_val          # already gzip bytes from pickle
            else:
                compressed = gzip.compress(model_val.encode())
            blob = numpy.frombuffer(compressed, dtype=numpy.uint8)
            dset.parent.create_dataset(
                family.accession + ".model",
                data=blob,
            )

        # Create links
        fam_link = f"/{group_path}/{family.accession}"
        if family.name:
            self.file.require_group(GROUP_LOOKUP_BYNAME)[str(family.name)] = (
                h5py.SoftLink(fam_link)
            )
        # In FamDB format version 0.5 we removed the /Families/ByAccession group as it's redundant
        # (all the data is in Families/<datasets> *and* HDF5 suffers from poor performance when
        # the number of entries in a group exceeds 200-500k.

        def add_stage_link(stage, accession):
            stage_group = self.file.require_group(GROUP_LOOKUP_BYSTAGE).require_group(
                stage.strip()
            )
            if accession not in stage_group:
                stage_group[accession] = h5py.SoftLink(fam_link)

        if family.search_stages:
            for stage in family.search_stages.split(","):
                add_stage_link(stage, family.accession)

        if family.buffer_stages:
            for stage in family.buffer_stages.split(","):
                stage = stage.split("[")[0]
                add_stage_link(stage, family.accession)

        LOGGER.debug(f"Added family {family.name} ({family.accession})")

    # Data Access Methods ------------------------------------------------------------------------------------------------
    def filter_stages(self, accession, stages):
        """Returns True if the family belongs to a search or buffer stage in 'stages'."""
        for stage in stages:
            grp = self.file[GROUP_LOOKUP_BYSTAGE].get(stage)
            if grp and accession in grp:
                return True

        return False

    # Family Getters --------------------------------------------------------------------------
    def get_family_by_accession(self, accession):
        """Returns the family with the given accession."""
        path = accession_bin(accession)
        if path in self.file:
            entry = self.file[path].get(accession)
            return get_family(entry)
        return None

    def get_family_by_name(self, name):
        """Returns the family with the given name."""
        # TODO: This will also suffer the performance issues seen with
        #       other groups that exceed 200-500k entries in a single group
        #       at some point.  This needs to be refactored to scale appropriately.
        # There are 24,768 names as of Dfam 3.8 - Anthony
        entry = self.file[GROUP_LOOKUP_BYNAME].get(name)
        if entry is None:
            name_lower = name.lower()
            for key in self.file[GROUP_LOOKUP_BYNAME].keys():
                if key.lower() == name_lower:
                    entry = self.file[GROUP_LOOKUP_BYNAME][key]
                    break
        return get_family(entry)


class FamDBRoot(FamDBLeaf):
    def __init__(self, filename, mode="r"):
        super(FamDBRoot, self).__init__(filename, mode)

        if mode == "r" or mode == "r+":
            names_ds = self.file.get(DATA_NAMES_CACHE)
            self.names_dump = json.loads(names_ds[()].decode()) if names_ds is not None else {}
            cache_ds = self.file.get(DATA_PARTITION_CACHE)
            self.partition_cache = json.loads(cache_ds[()].decode()) if cache_ds is not None else {}
            self.file_info = self.get_file_info()
            self.__lineage_cache = {}

    @FamDBLeaf._change_logger
    def write_full_taxonomy(self, tax_db):
        """
        Takes a map of TaxaNodes keyed by tax_id.
        Writes taxonomy nodes to the database including parent-child relationships.
        Also caches all taxa names as a node:[names] JSON string loaded at __init__.
        Partition assignments are no longer stored per-node; use write_partition_cache()
        after all component files are built.
        """
        LOGGER.debug(f"Writing Full Taxonomy Tree Root File")
        start = time.perf_counter()

        names_dump = {}
        count = 0
        for node in tax_db:
            count += 1
            group = self.file.require_group(GROUP_NODES).require_group(
                str(tax_db[node].tax_id)
            )
            parent_id = int(tax_db[node].parent_id) if node != 1 else None
            if parent_id:
                group.create_dataset(DATA_PARENT, data=numpy.array([parent_id]))

            child_ids = []
            for child in tax_db[node].children:
                child_ids += [int(child.tax_id)]
            group.create_dataset(DATA_CHILDREN, data=numpy.array(child_ids))

            names = tax_db[node].names
            group.create_dataset(DATA_TAXANAMES, data=numpy.array(names, dtype="S"))
            names_dump[node] = names

        LOGGER.debug(f"Writing Name Cache String")
        self.file.create_dataset(
            DATA_NAMES_CACHE, data=numpy.array(json.dumps(names_dump), dtype="S")
        )

        delta = time.perf_counter() - start
        LOGGER.info(f"Wrote {count} taxonomy nodes in full tree in {delta:.1f}s")

    def update_pruned_taxa(self, tree):
        """
        Takes a map of TaxaNodes
        Updates the nodes to include sparse parent-child relationships
        based on which nodes have family data associated with them
        """
        for id in tree:
            node = tree[id]
            val_children = [int(child) for child in node.val_children]
            val_parent = int(node.val_parent) if node.val_parent else None
            group = self.file[GROUP_NODES][str(id)]
            if group.get(DATA_VAL_CHILDREN):
                del group[DATA_VAL_CHILDREN]
            group.create_dataset(
                DATA_VAL_CHILDREN,
                data=numpy.array(val_children),
                shape=(len(val_children),),
                dtype="i8",
            )
            if val_parent:
                if group.get(DATA_VAL_PARENT):
                    del group[DATA_VAL_PARENT]
                group.create_dataset(
                    DATA_VAL_PARENT,
                    data=numpy.array([val_parent]),
                    shape=(1,),
                    dtype="i8",
                )

    @FamDBLeaf._change_logger
    def write_repeatpeps(self, infile):
        """
        Writing RepeatPeps to its own group as one big string.
        For now, only RepeatModeler consumes this, and does so
        by loading the whole file, so no need to do more
        """
        LOGGER.info(f"Writing RepeatPeps From File: {infile}")
        fasta = is_fasta(infile)
        if fasta:
            with open(infile, "r") as file:
                repeatpeps_str = file.read()
                rp_data = self.file.create_dataset(
                    GROUP_REPEATPEPS, shape=1, dtype=h5py.string_dtype()
                )
                rp_data[:] = repeatpeps_str
            LOGGER.info("RepeatPeps Saved")
        else:
            LOGGER.error(f"File {infile} not in FASTA format, write cancelled")

    def get_repeatpeps(self):
        """
        Retrieve RepeatPeps File
        """
        return self.file.get(GROUP_REPEATPEPS)[0].decode(
            encoding="UTF-8", errors="strict"
        )

    # currently unused:
    # def get_family_names(self):
    #     """Returns a list of names of families in the database."""
    #     return sorted(self.file[GROUP_LOOKUP_BYNAME].keys(), key=str.lower)

    def get_taxon_names(self, tax_id):
        """
        Checks names_dump for each partition and returns a list of [name_class, name_value, partition]
        of the taxon given by 'tax_id'.
        """
        nodes = self.file[GROUP_NODES]
        node = nodes.get(str(tax_id))
        if node:
            return [
                [name.decode() for name in name_pair]
                for name_pair in node[DATA_TAXANAMES][:]
            ]
        return []

    def get_taxon_name(self, tax_id, kind="scientific name"):
        """
        Returns the first name of the given 'kind' for the taxon given by 'tax_id',
        along with the component partition dict for that taxon.
        Returns ("Not Found", "N/A") if the taxon is not found.
        """
        failure = ("Not Found", "N/A")

        nodes = self.file[GROUP_NODES]
        node = nodes.get(str(tax_id))
        if not node:
            return failure

        names = [
            [name.decode() for name in name_pair]
            for name_pair in node[DATA_TAXANAMES][:]
        ]
        partition = self.partition_cache.get(str(tax_id))

        if names:
            for name in names:
                if name[0] == kind:
                    return [name[1], partition]
        return failure

    def search_taxon_names(self, text, kind=None, search_similar=False):
        """
        Searches 'self' for taxons with a name containing 'text', returning an
        iterator that yields a tuple of (id, is_exact, partition) for each matching node.
        Each id is returned at most once, and if any of its names are an exact
        match the whole node is treated as an exact match.

        If 'similar' is True, names that sound similar will also be considered
        eligible.

        A list of strings may be passed as 'kind' to restrict what kinds of
        names will be searched.
        """

        text = text.lower()
        for tax_id, names in self.names_dump.items():
            matches = False
            exact = False
            for name_cls, name_txt in names:
                name_txt = name_txt.lower()
                if kind is None or kind == name_cls:
                    if text == name_txt:
                        matches = True
                        exact = True
                    elif name_txt.startswith(text + " <"):
                        matches = True
                        exact = True
                    elif text == sanitize_name(name_txt):
                        matches = True
                        exact = True
                    elif text in name_txt:
                        matches = True
                    elif search_similar and sounds_like(text, name_txt):
                        matches = True

            if matches:
                partition = self.find_taxon(tax_id)
                yield [int(tax_id), exact, partition]

    def resolve_species(self, term, kind=None, search_similar=False):
        """
        Resolves 'term' as a species or clade in 'self'. If 'term' is a number,
        it is a taxon id. Otherwise, it will be searched for in 'self' in the
        name fields of all taxa. A list of strings may be passed as 'kind' to
        restrict what kinds of names will be searched.

        If 'search_similar' is True, a "sounds like" search will be tried
        first. If it is False, a "sounds like" search will still be performed

        if no results were found.

        This function returns a list of tuples (taxon_id, is_exact) that match
        the query. The list will be empty if no matches were found.
        """
        # Try as a number
        try:
            tax_id = int(term)
            if str(tax_id) in self.names_dump:
                partition = self.find_taxon(tax_id)
                return [[tax_id, partition, True]]

            return []
        except ValueError:
            pass

        # Perform a search by name, splitting between exact and inexact matches for sorting
        exact = []
        inexact = []
        for tax_id, is_exact, partition in self.search_taxon_names(
            term, kind, search_similar
        ):
            hit = [tax_id, partition]
            if is_exact:
                exact += [hit]
            else:
                inexact += [hit]

        # Combine back into one list, with exact matches first
        results = [[*hit, True] for hit in exact]
        for hit in inexact:
            results += [[*hit, False]]

        if len(results) == 0 and not search_similar:
            # Try a sounds-like search (currently soundex)
            similar_results = self.resolve_species(term, kind, True)
            if similar_results:
                print(
                    "No results were found for that name, but some names sound similar:",
                    file=sys.stderr,
                )
                for tax_id, partition, exact in similar_results:
                    names = self.get_taxon_names(tax_id)
                    print(
                        tax_id,
                        ", ".join([f"{n}" for n in names]),
                        file=sys.stderr,
                    )

        return results

    def resolve_one_species(self, term, kind=None):
        """
        Resolves 'term' in 'dbfile' as a taxon id or search term unambiguously.
        Parameters are as in the 'resolve_species' method.
        Returns None if not exactly one result is found,
        and prints details to the screen.
        """

        results = self.resolve_species(term, kind)

        # Check for a single exact match first, to any field
        exact_matches = []

        for result in results:  # result -> [tax_id, partition, exact]
            if result[2]:
                exact_matches += [[result[0], result[1]]]
        if len(exact_matches) == 1:
            return exact_matches[0]

        if len(results) == 1:
            return results[0][:2]
        elif len(results) > 1:
            print(
                f"""Ambiguous search term '{term}' (found {len(results)} results, {len(exact_matches)} exact).
Please use a more specific name or taxa ID, which can be looked
up with the 'names' command.""",
                file=sys.stderr,
            )
            return "Ambiguous", "Ambiguous"
        return None, None

    def get_sanitized_name(self, tax_id):
        """
        Returns the "sanitized name" of tax_id, which is the sanitized version
        of the scientific name.
        Used in EMBL exports
        """

        name = self.get_taxon_name(tax_id, "scientific name")
        if name:
            name = sanitize_name(name[0])
        return name

    def get_lineage(self, tax_id, **kwargs):
        """
        Returns the lineage of 'tax_id'. Recognized kwargs: 'descendants' to include
        descendant taxa, 'ancestors' to include ancestor taxa.
        IDs are returned as a nested list, for example
        [ 1, [ 2, [3, [4]], [5], [6, [7]] ] ]
        where '2' may have been the passed-in 'tax_id'.
        """

        group_nodes = self.file[GROUP_NODES]
        ancestors = True if kwargs.get("ancestors") else False
        descendants = True if kwargs.get("descendants") else False
        children_key = (
            DATA_VAL_CHILDREN if not kwargs.get("complete") else DATA_CHILDREN
        )
        parent_key = DATA_VAL_PARENT if not kwargs.get("complete") else DATA_PARENT

        if descendants:

            def descendants_of(tax_id):
                descendants = [
                    int(tax_id)
                ]  # h5py is based on numpy, need to cast numpy base64 to python int for serialization in Lineage class
                for child in group_nodes[str(tax_id)][children_key]:
                    descendants += [descendants_of(child)]
                return descendants

            tree = descendants_of(tax_id)
        else:
            tree = [tax_id]

        if ancestors:
            while tax_id:
                node = group_nodes[str(tax_id)]
                if parent_key in node:
                    tax_id = node[parent_key][0]
                    tree = [
                        int(tax_id),
                        tree,
                    ]  # h5py is based on numpy, need to cast numpy base64 to python int for serialization in Lineage class
                else:
                    tax_id = None

        return tree

    def get_lineage_path(self, tax_id, cache=True, partition=True, complete=False):
        """
        Returns a list of strings encoding the lineage for 'tax_id'.
        """

        if cache and tax_id in self.__lineage_cache:
            return self.__lineage_cache[tax_id]
        tree = self.get_lineage(tax_id, ancestors=True, complete=complete)
        lineage = []

        while tree:
            node = tree[0]
            if len(tree) > 1:
                found = False
                for t in tree[1:]:
                    if type(t) == list:
                        tree = t
                        found = True
                        break
                if not found:
                    tree = None
            else:
                tree = None

            tax_name = self.get_taxon_name(node, "scientific name")
            if not partition:
                tax_name = tax_name[0]
            lineage += [tax_name]

        if cache:
            self.__lineage_cache[tax_id] = lineage

        return lineage

    def get_partition_for_taxon(self, tax_id, component):
        """
        Returns the partition number (int) for 'tax_id' in the given component,
        or None if the taxon has no families of that component type.
        'component' must be one of COMPONENT_CC, COMPONENT_CH, COMPONENT_UC, COMPONENT_UH.
        """
        entry = self.partition_cache.get(str(tax_id))
        if entry:
            return entry.get(component)
        return None

    def find_taxon(self, tax_id):
        """
        Returns a dict mapping component type to partition number for the given taxon.
        e.g. {"cc": 0, "ch": 2, "uc": 1, "uh": 1}
        Values are None for components with no families at this taxon.
        Returns None if taxon is not in the partition cache at all.
        """
        return self.partition_cache.get(str(tax_id))

    def get_families_for_taxon(self, tax_id, curated_only=False, uncurated_only=False):
        """
        Returns a list of the accessions for each family directly associated with 'tax_id'.
        Reads from the root file's Lookup/ByTaxon, which covers all families.
        """
        key = f"{GROUP_LOOKUP_BYTAXON}/{tax_id}"
        if key not in self.file:
            return []
        taxon_group = self.file[key]
        # New format: group contains an 'accessions' dataset
        if "accessions" in taxon_group:
            accessions = list(taxon_group["accessions"].asstr()[()])
        else:
            # Legacy group-of-subgroups format
            accessions = list(taxon_group.keys())

        if curated_only:
            return list(filter(lambda a: filter_curated(a, True), accessions))
        elif uncurated_only:
            return list(filter(lambda a: filter_curated(a, False), accessions))
        else:
            return accessions

    @FamDBLeaf._change_logger
    def write_lookup_bytaxon(self, family_taxon_map):
        """
        Writes Lookup/ByTaxon to the root file.
        family_taxon_map: {tax_id_str: [accession, ...]}
        Called as a post-processing step after all component files are built.
        All families (curated + uncurated) are included.

        Accessions are stored as a variable-length string dataset under each
        taxon group - one dataset write per taxon instead of one sub-group
        per accession.  get_families_for_taxon() reads the dataset back.
        """
        LOGGER.info("Writing Lookup/ByTaxon to root file")
        start = time.perf_counter()
        count = 0
        str_dtype = h5py.string_dtype()
        bytaxon_group = self.file.require_group(GROUP_LOOKUP_BYTAXON)
        for tax_id, accessions in family_taxon_map.items():
            taxon_group = bytaxon_group.require_group(str(tax_id))
            taxon_group.create_dataset(
                "accessions",
                data=numpy.array(accessions, dtype=object),
                dtype=str_dtype,
            )
            count += len(accessions)
        delta = time.perf_counter() - start
        LOGGER.info(f"Wrote {count} taxon-family entries in {delta:.2f}s")

    @FamDBLeaf._change_logger
    def write_partition_cache(self, partition_cache):
        """
        Writes the PartitionCache JSON blob to the root file.
        partition_cache: {tax_id_str: {"cc": int|None, "ch": int|None,
                                        "uc": int|None, "uh": int|None}}
        Only taxa with at least one non-None component entry are included.
        """
        LOGGER.info("Writing PartitionCache to root file")
        if self.file.get(DATA_PARTITION_CACHE):
            del self.file[DATA_PARTITION_CACHE]
        self.file.create_dataset(
            DATA_PARTITION_CACHE,
            data=numpy.array(json.dumps(partition_cache), dtype="S"),
        )
        LOGGER.info(f"PartitionCache written ({len(partition_cache)} taxa)")

    def _add_family_taxon_links(self, accession, clade_ids):
        """
        Add a single family->taxon mapping to Lookup/ByTaxon.
        Called during append to update the root index without a full rewrite.
        No changelog entry is created (individual additions are logged at the
        FamDB level).
        """
        str_dtype = h5py.string_dtype()
        bytaxon_group = self.file.require_group(GROUP_LOOKUP_BYTAXON)
        for clade_id in clade_ids:
            taxon_group = bytaxon_group.require_group(str(clade_id))
            if "accessions" in taxon_group:
                existing = list(taxon_group["accessions"].asstr()[()])
                if accession not in existing:
                    existing.append(accession)
                    del taxon_group["accessions"]
                    taxon_group.create_dataset(
                        "accessions",
                        data=numpy.array(existing, dtype=object),
                        dtype=str_dtype,
                    )
            else:
                taxon_group.create_dataset(
                    "accessions",
                    data=numpy.array([accession], dtype=object),
                    dtype=str_dtype,
                )

    def get_all_taxa_names(self):
        """
        Returns all taxa names in database.
        Names are cached as taxa : [[name type, name], ...]
        Names are returned as {sanitized_lowercase_name: taxa...}
        Used for mapping EMBL file names to taxa nodes
        Used in append command
        """
        return {
            name[1].lower(): taxon
            for taxon, names in self.names_dump.items()
            for name in names
            if name[0] == "sanitized scientific name" or name[0] == "sanitized synonym"
        }


class FamDB:

    def __init__(self, db_dir, mode, exclude=[]):
        """
        Initialize from a directory containing a v3 famdb dataset.

        File naming convention:
          <base>.0.h5                              - root file (taxonomy index)
          <base>.<curated|uncurated>.<consensus|hmm>.<N>.h5  - component files

        self.files[0]        = FamDBRoot (root file)
        self.components      = {component_type: {partition_num: FamDBLeaf}}
        self.files           = {0: root} union {all component leaf files keyed by a
                                 unique int id} for backward-compat with internal
                                 helpers that iterate self.files.
        """
        self.files = {}
        # {component_type: {partition_num: FamDBLeaf}}
        self.components = {ct: {} for ct in COMPONENT_TYPES}

        h5_files = sorted(os.listdir(db_dir))
        root_file = None
        db_prefix = None

        # First pass: locate the root file and determine the base prefix
        for filename in h5_files:
            # Skip component files - they also match the root regex when partition=0
            if FAMDB_COMPONENT_FILE_RE.match(filename):
                continue
            m = FAMDB_ROOT_FILE_RE.match(filename)
            if m:
                if db_prefix is not None and m.group(1) != db_prefix:
                    LOGGER.error(
                        "Multiple famdb root files found in " + db_dir +
                        ". Each famdb database should be in a separate folder."
                    )
                    exit(1)
                db_prefix = m.group(1)
                root_file = filename

        if root_file is None:
            if h5_files:
                LOGGER.error(
                    "A famdb root file (*.0.h5) is not present in " + db_dir
                )
            else:
                LOGGER.error("No .h5 files found in " + db_dir)
            exit(1)

        self.files[0] = FamDBRoot(f"{db_dir}/{root_file}", mode)

        # Second pass: load component files
        _file_id_counter = 1  # unique ints for self.files backward compat
        for filename in h5_files:
            m = FAMDB_COMPONENT_FILE_RE.match(filename)
            if not m:
                continue
            prefix, curated_str, model_str, part_str = m.groups()
            if prefix != db_prefix:
                continue
            part_num = int(part_str)
            # Map (curated|uncurated, consensus|hmm) -> component type
            component_type = (
                COMPONENT_CC if curated_str == "curated" and model_str == "consensus" else
                COMPONENT_CH if curated_str == "curated" and model_str == "hmm" else
                COMPONENT_UC if curated_str == "uncurated" and model_str == "consensus" else
                COMPONENT_UH
            )
            if component_type in exclude:
                continue
            leaf = FamDBLeaf(f"{db_dir}/{filename}", mode)
            self.components[component_type][part_num] = leaf
            self.files[_file_id_counter] = leaf
            _file_id_counter += 1

        file_info = self.files[0].get_file_info()
        self.db_dir = db_dir
        self.file_map = file_info[META_FILE_MAP]
        self.uuid = file_info[META_META][META_UUID]
        self.db_version = file_info[META_META][META_DB_VERSION]
        self.db_date = file_info[META_META][META_DB_DATE]

        # Validate all component files match root metadata
        partition_err_files = []
        for fid, fobj in self.files.items():
            if fid == 0:
                continue
            meta = fobj.get_file_info()[META_META]
            if (
                self.uuid != meta[META_UUID]
                or self.db_version != meta[META_DB_VERSION]
                or self.db_date != meta[META_DB_DATE]
            ):
                partition_err_files.append(fobj.filename)
        if partition_err_files:
            LOGGER.error(f"Files From Different Partitioning Runs: {partition_err_files}")
            exit()

        change_err_files = []
        for fobj in self.files.values():
            if fobj.interrupt_check():
                change_err_files.append(fobj.filename)
        if change_err_files:
            LOGGER.error(f"Files Interrupted During Edit: {change_err_files}")
            exit()

    def _check_component(self, component_type):
        """
        Returns True if at least one file of the given component type is loaded.
        Logs a warning if the component is missing, describing what to download.
        """
        if self.components.get(component_type):
            return True
        LOGGER.warning(
            f"Component '{component_type}' is not installed in {self.db_dir}. "
            f"The requested query requires this data. "
            f"Please download the corresponding files from Dfam."
        )
        return False

    # Data writing methods ---------------------------------------------------------------------------------------
    def _make_tax_node(self, id, hdf5_node=None, value=False, read_pruned=False):
        """
        Build a TaxNode from HDF5 node data.

        id         : taxonomy id (string or int)
        hdf5_node  : already-fetched HDF5 group; looked up from GROUP_NODES if None
        value      : whether this node has associated family data (val flag)
        read_pruned: if True, also populate val_children and val_parent from HDF5
        """
        if hdf5_node is None:
            hdf5_node = self.files[0].file[GROUP_NODES][str(id)]
        children = hdf5_node[DATA_CHILDREN][()] if hdf5_node[DATA_CHILDREN].size > 0 else []
        parent = (
            hdf5_node[DATA_PARENT][()][0]
            if hdf5_node.get(DATA_PARENT) and hdf5_node[DATA_PARENT].size > 0
            else None
        )
        node = TaxNode(id, str(parent) if parent else None)
        node.val = value
        node.children = children
        if read_pruned:
            node.val_children = (
                hdf5_node[DATA_VAL_CHILDREN][()] if hdf5_node[DATA_VAL_CHILDREN].size > 0 else []
            )
            node.val_parent = (
                hdf5_node[DATA_VAL_PARENT][()][0]
                if hdf5_node.get(DATA_VAL_PARENT) and hdf5_node[DATA_VAL_PARENT].size > 0
                else None
            )
        return node

    def build_pruned_tree(self):
        """
        Establishes a sparse taxonomy tree where parent-child relationships are restricted to
        nodes with associated family data. For example, a node will be assigned a sparse parent
        as the closest ancestor node with data, rather than its actual parent node, if its
        actual parent node is empty.
        This method exists in FamDB instead of FamDBRoot because it is subject to change after an append,
        and because the associated data is stored in FamDBLeaf files
        """

        def traverse_val_parents(tree, id):
            """Recurse up the tree ancestor by ancestor until it finds the nearest ancestor with data"""
            node = tree[id]
            if node.parent_id:
                parent = tree.get(node.parent_id)
                if parent:
                    if parent.val:
                        return parent.tax_id
                    else:
                        return traverse_val_parents(tree, parent.tax_id)
            else:
                return None

        def traverse_val_children(tree, id, node_id):
            """
            Adds node to it's parent's list of sparse children
            Continues recursion until it finds an ancestor with data
            """
            node = tree[id]
            if node.parent_id:
                parent = tree.get(node.parent_id)
                if parent:
                    parent.val_children += [node_id]
                    if not parent.val:
                        traverse_val_children(tree, parent.tax_id, node_id)

        LOGGER.info("Reading Taxonomy Tree")
        # read taxonomy tree
        tree = {
            node: self.files[0].file[GROUP_NODES][node]
            for node in self.files[0].file[GROUP_NODES]
        }
        LOGGER.info("Determining Which Nodes Have Associated Families")
        # In v3 Lookup/ByTaxon lives in the root file only and covers all families
        bytaxon = self.files[0].file[GROUP_LOOKUP_BYTAXON]
        vals = set(
            tax_id
            for tax_id in bytaxon
            if bool(bytaxon[tax_id].keys())
        )

        # build TaxNodes in tree
        for id in tree:
            tree[id] = self._make_tax_node(id, hdf5_node=tree[id], value=(id in vals))

        LOGGER.info("Full Tree Prepared")
        # assign each node a val_parent
        for id in tree:
            node = tree[id]
            node.val_parent = traverse_val_parents(tree, node.tax_id)

        # add each node with a value to it's parents as a val_child
        for id in tree:
            node = tree[id]
            if node.val:
                traverse_val_children(tree, node.tax_id, node.tax_id)

        LOGGER.info("Pruned Tree Prepared")

        # update database nodes - changelog only on root, not all leaf files
        message = "Pruned Tree Written"
        ts = self.files[0].update_changelog(message)
        self.files[0].update_pruned_taxa(tree)
        self.files[0]._verify_change(ts, message)
        LOGGER.info(message)

    def rebuild_pruned_tree(self, new_val_taxa):
        """
        This method takes a list/set of taxon ids that did not have families associated with them,
        but do now due to a recent append command. It resets the val_parent/val_child links in the
        taxonomy tree to account for the fact that there is new data in the tree.
        It assumes that a subject node's val_parent and all ancestor nodes between them will set it
        as one of thier val_children in place of any val_children that it used to share with the
        subject node.
        Likewise, it assumes that any of it's val_children and all child nodes betweem will replace
        thier val_parents with the subject node.
        """

        def build_taxa_node(id, value=False):
            return self._make_tax_node(id, value=value, read_pruned=True)

        # RMH: This parameter default pattern "foo=[]" is dangerous.  The
        #      list generated is global and gets reused between independent
        #      invocations!
        # def climb_non_val_parents(node, ancestor_path=[]):
        #    """collects the nodes between a node and it's val_parent, not inclusive"""
        #    if node.parent_id != node.val_parent:
        #        parent_node = build_taxa_node(node.parent_id)
        #        ancestor_path += [parent_node]
        #        climb_non_val_parents(parent_node, ancestor_path)
        #    return ancestor_path

        def climb_non_val_parents(target_id, node, ancestor_path=None):
            """Collects TaxNodes between a given node and a ancestral
            node defined by target_id (exclusive)."""

            if ancestor_path == None:
                ancestor_path = []

            if hasattr(node, "parent_id") and str(node.parent_id) != str(target_id):
                parent_node = build_taxa_node(node.parent_id)
                ancestor_path += [parent_node]
                ancestor_path = climb_non_val_parents(
                    target_id, parent_node, ancestor_path
                )
            return ancestor_path

        message = "Pruned Tree Updated"
        rec = self.append_start_changelog(message)
        update_nodes = {}
        for id in new_val_taxa:
            node = build_taxa_node(id, value=True)

            seen_tax_ids = set()
            # Fix the descendants of the target node
            for val_child in node.val_children:
                child_node = build_taxa_node(val_child, value=True)

                # Fix the val_child node itself
                if child_node.tax_id not in seen_tax_ids:
                    seen_tax_ids.add(child_node.tax_id)
                    child_node.val_parent = id
                    update_nodes[child_node.tax_id] = child_node

                # Fix the ancestors up until (but not including) the target node
                for ancestor in climb_non_val_parents(id, child_node):
                    if ancestor.tax_id not in seen_tax_ids:
                        seen_tax_ids.add(ancestor.tax_id)
                        ancestor.val_parent = id
                        update_nodes[ancestor.tax_id] = ancestor

            # Gather all nodes above the target node up until its val_parent
            if node.val_parent:
                change_ancestors = [build_taxa_node(node.val_parent, value=True)]
                change_ancestors += climb_non_val_parents(node.val_parent, node)

                # All nodes above it should point to it as well, instead of any of its val_children
                for ansc_node in change_ancestors:
                    # remove any val_children that are below this node
                    for vid in node.val_children:
                        ansc_node.val_children = ansc_node.val_children[
                            ansc_node.val_children != vid
                        ]
                    # add this node to the ancestral val_children
                    ansc_node.val_children = numpy.append(ansc_node.val_children, id)
                    update_nodes[ansc_node.tax_id] = ansc_node

            # update the tree for each newly val'd taxon, to avoid tangling pointers when multiple updates occur on the same path
            self.files[0].update_pruned_taxa(update_nodes)
            update_nodes = {}
        self.append_finish_changelog(message, rec)
        LOGGER.info(message)

    def set_db_info(self, name, version, date, desc, copyright_text):
        """Method for resetting metadata"""
        for fobj in self.files.values():
            partition_num = fobj.get_partition_num()
            file_info = fobj.get_file_info()
            fobj.set_metadata(partition_num, file_info, name, version, date, copyright_text)
            fobj.update_description(desc)

    def append_start_changelog(self, message):
        """Called when an append command starts"""
        rec = {}
        for fid, fobj in self.files.items():
            rec[fid] = fobj.update_changelog(message)
        return rec

    def append_finish_changelog(self, message, rec):
        """Called when an append command finishes successfully"""
        for fid, timestamp in rec.items():
            self.files[fid]._verify_change(timestamp, message)

    def update_changelog(self, added_ctr, total_ctr, file_counts, infile):
        """Used to add a context log after an append command"""
        filename = infile.split("/")[-1]
        for fid, fobj in self.files.items():
            if fid in file_counts:
                fobj.update_changelog(
                    f"Added {file_counts[fid]} of {total_ctr} Families From {filename}",
                    verified=True,
                )
            else:
                fobj.update_changelog(
                    f"Found No Relevant Families From {filename}", verified=True
                )
            if fid == 0:
                fobj.update_changelog(
                    f"Total Families {added_ctr} of {total_ctr} Added To Local Files From {filename}",
                    verified=True,
                )

    # Data access methods ---------------------------------------------------------------------------------------
    def show_files(self):
        """Shows loaded component files and their counts."""
        print(f"\nInstalled Components\n--------------------")
        component_labels = {
            COMPONENT_CC: "Curated Consensus",
            COMPONENT_CH: "Curated HMMs",
            COMPONENT_UC: "Uncurated Consensus",
            COMPONENT_UH: "Uncurated HMMs",
        }
        ct_key = {COMPONENT_CC: "cc", COMPONENT_CH: "ch",
                  COMPONENT_UC: "uc", COMPONENT_UH: "uh"}
        for ct in COMPONENT_TYPES:
            label = component_labels[ct]
            print(f"\n {label}:")
            # All partitions expected for this component type, from the file map
            expected = {
                int(k.split(".")[1]): v
                for k, v in self.file_map.items()
                if k.startswith(ct_key[ct] + ".")
            }
            if not expected:
                print(f"     [ Not Installed ]")
                continue
            is_hmm = ct in (COMPONENT_CH, COMPONENT_UH)
            rows = []
            for part_num in sorted(expected):
                fm_entry = expected[part_num]
                root_name = fm_entry.get("T_root_name", "")
                leaf = self.components[ct].get(part_num)
                if leaf is not None:
                    counts = leaf.get_counts()
                    count = counts["hmm"] if is_hmm else counts["consensus"]
                    rows.append((part_num, fm_entry["filename"], root_name, count))
                else:
                    rows.append((part_num, fm_entry["filename"], root_name, None))
            max_prefix = max(len(f"     partition {r[0]} [{r[1]}]:") for r in rows)
            max_name = max(len(r[2]) for r in rows)
            present = [r[3] for r in rows if r[3] is not None]
            max_num = max((len(f"{c:,}") for c in present), default=3)
            for part_num, filename, root_name, count in rows:
                prefix = f"     partition {part_num} [{filename}]:"
                if count is not None:
                    print(f"{prefix:<{max_prefix}}  {root_name:<{max_name}}  {count:>{max_num},} families")
                else:
                    print(f"{prefix:<{max_prefix}}  {root_name:<{max_name}}  --- not present ---")
        print()

    def show_history(self):
        """Iterates over all present files and prints each history"""
        print(f"\nFile History\n-----------------")
        for fobj in self.files.values():
            print(fobj.get_history())

    def print_info(self, history=False):
        """Print stored metadata, file inventory, and optionally changelog history."""
        db_info = self.get_metadata()
        counts = self.get_counts()
        print()
        print(
            f"""\
FamDB Directory               : {os.path.realpath(self.db_dir)}
FamDB Creation Format Version : {db_info["famdb_version"]}
FamDB Creation Date           : {db_info["created"]}

Database : {db_info["name"]}
Version  : {db_info["db_version"]}
Date     : {db_info["date"]}

{db_info["description"]}

{counts['file']} Partitions Present
Total consensus sequences present: {counts["consensus"]}
Total HMMs present               : {counts["hmm"]}
"""
        )
        self.show_files()
        if history:
            self.show_history()

    def get_counts(self):
        """Method gets collected counts from each file present"""
        counts = {"consensus": 0, "hmm": 0, "file": 0}
        for fobj in self.files.values():
            file_counts = fobj.get_counts()
            counts["consensus"] += file_counts["consensus"]
            counts["hmm"] += file_counts["hmm"]
            counts["file"] += 1
        return counts

    def assemble_filters(self, **kwargs):
        """Define family filters (logically ANDed together)"""
        filters = []
        if kwargs.get("curated_only"):
            filters += [lambda a, f: filter_curated(a, True)]
        if kwargs.get("uncurated_only"):
            filters += [lambda a, f: filter_curated(a, False)]

        filter_stage = kwargs.get("stage")
        stages = []
        if filter_stage is not None:
            if filter_stage == 0:
                # RMH: 6/27/25
                # stage 0 = 'no stage defined'
                # so filter out anything with a defined search stage
                filters += [lambda a, f: filter_defined_search_stages(f())]
            elif filter_stage == 80:
                # "stage 80" = "all stages", so skip filtering
                pass
            elif filter_stage == 95:
                # "stage 95" = this specific stage list:
                stages = ["35", "50", "55", "60", "65", "70", "75"]
                filters += [lambda a, f: self.filter_stages(a, stages)]
            else:
                stages = [str(filter_stage)]
                filters += [lambda a, f: self.filter_stages(a, stages)]

        # HMM only: add a search stage filter to "un-list" families that were
        # allowed through only because they match in buffer stage
        if kwargs.get("is_hmm") and stages:
            filters += [lambda a, f: filter_search_stages(f(), stages)]

        repeat_type = kwargs.get("repeat_type")
        if repeat_type:
            repeat_type = repeat_type.lower()
            filters += [lambda a, f: filter_repeat_type(f(), repeat_type)]

        name = kwargs.get("name")
        if name:
            name = name.lower()
            filters += [lambda a, f: filter_name(f(), name)]

        return filters, stages, repeat_type, name

    def get_accessions_filtered(self, **kwargs):
        """
        Returns an iterator that yields accessions for the given search terms.

        Filters are specified in kwargs:
            tax_id: int
            ancestors: boolean, default False
            descendants: boolean, default False
                If none of (tax_id, ancestors, descendants) are
                specified, *all* families will be checked.
            curated_only = boolean
            uncurated_only = boolean
            stage = int
            is_hmm = boolean
            repeat_type = string (prefix)
            name = string (prefix)
                If any of stage, repeat_type, or name are
                omitted (or None), they will not be used to filter.

                The behavior of 'stage' depends on 'is_hmm': if is_hmm is True,
                stage must match in SearchStages (a match in BufferStages is not
                enough).
        """

        if not ("tax_id" in kwargs or "ancestors" in kwargs or "descendants" in kwargs):
            tax_id = 1
            ancestors = True
            descendants = True
        else:
            tax_id = kwargs["tax_id"]
            ancestors = kwargs.get("ancestors") or False
            descendants = kwargs.get("descendants") or False

        filters, stages, repeat_type, name_filter = self.assemble_filters(**kwargs)

        # Recursive iterator flattener
        def walk_tree(tree):
            """Returns all elements in 'tree' with all levels flattened."""
            if hasattr(tree, "__iter__"):
                for elem in tree:
                    yield from walk_tree(elem)
            else:
                yield tree

        seen = set()

        curated_only = kwargs.get("curated_only", False)
        uncurated_only = kwargs.get("uncurated_only", False)

        def iterate_accs():
            # special case: Searching the whole database in a specific
            # stage only is a common usage pattern in RepeatMasker.
            # When searching the whole database instead of a species,
            # the number of accessions to read through is shorter
            # when going off of only the stage indexes.
            files = self.files
            if (
                tax_id == 1
                and descendants
                and stages
                and not repeat_type
                and not name_filter
            ):
                for stage in stages:
                    for file in files:
                        by_stage = files[file].file.get(GROUP_LOOKUP_BYSTAGE)
                        if by_stage:
                            grp = by_stage.get(stage)
                            if grp:
                                yield from grp.keys()

            # special case: Searching the whole database, going directly via
            # Families/ is faster than repeatedly traversing the tree
            elif tax_id == 1 and descendants:
                for file in files:
                    if GROUP_FAMILIES not in files[file].file:
                        continue
                    names = families_iterator(
                        files[file].file[GROUP_FAMILIES], GROUP_FAMILIES
                    )
                    for name in names:
                        yield name
            else:
                lineage = self.get_lineage(
                    tax_id, ancestors=ancestors, descendants=descendants
                )
                for node in walk_tree(lineage):
                    fams = self.get_families_for_taxon(
                        node, curated_only=curated_only, uncurated_only=uncurated_only
                    )
                    if fams:
                        yield from fams

        for accession in iterate_accs():
            if accession in seen:
                continue
            seen.add(accession)

            cached_family = None

            def family_getter():
                nonlocal cached_family
                if not cached_family:
                    path = accession_bin(accession)
                    for file in self.files:
                        if self.files[file].file.get(path):
                            fam = self.files[file].file[path].get(accession)
                            if fam:
                                cached_family = fam
                return cached_family

            match = True
            for filt in filters:
                if not filt(accession, family_getter):
                    match = False
            if match:
                yield accession

    def fasta_all(self, group):
        """
        Collects all families in a Families sub-group (e.g. '/DF', '/Aux').
        Searches curated and uncurated consensus component files only, since
        consensus sequences are what fasta_all callers need.
        """
        seen = set()
        consensus_files = list(self.components[COMPONENT_CC].values()) + \
                          list(self.components[COMPONENT_UC].values())
        for leaf in consensus_files:
            target = GROUP_FAMILIES + group
            if target in leaf.file:
                for name in families_iterator(leaf.file[target], target):
                    if name not in seen:
                        seen.add(name)
                        yield self.get_family_by_accession(name)

    # Root Wrapper methods ---------------------------------------------------------------------------------------
    def resolve_names(self, term):
        """Method to find names matching the search term and map them to the correct file"""
        entries = []
        for tax_id, partition, is_exact in self.files[0].resolve_species(term):
            names = self.files[0].get_taxon_names(tax_id)
            entries += [[tax_id, is_exact, partition, names]]
        return entries

    def get_lineage_path(self, tax_id, **kwargs):
        """method used in EMBL exports"""
        partition = (
            kwargs.get("partition") if kwargs.get("partition") is not None else True
        )
        cache = kwargs.get("cache") if kwargs.get("cache") is not None else True
        complete = (
            kwargs.get("complete") if kwargs.get("complete") is not None else True
        )
        return self.files[0].get_lineage_path(
            tax_id, cache=cache, partition=partition, complete=complete
        )

    def get_sanitized_name(self, tax_id):
        """Wrapper method for the Root get_sanitized_name method"""
        return self.files[0].get_sanitized_name(tax_id)

    def get_lineage(self, tax_id, **kwargs):
        """Wrapper method for the Root get_lineage method"""
        return self.files[0].get_lineage(tax_id, **kwargs)

    def resolve_one_species(self, term):
        """Wrapper method for the Root resolve_one_species method"""
        return self.files[0].resolve_one_species(term)

    def get_metadata(self):
        """Wrapper method for the Root get_metadata method"""
        return self.files[0].get_metadata()

    def get_taxon_name(self, tax_id, kind):
        """Wrapper method for the Root get_taxon_name method"""
        return self.files[0].get_taxon_name(tax_id, kind)

    def find_taxon(self, tax_id):
        """Wrapper method for the Root find_taxon method"""
        return self.files[0].find_taxon(tax_id)

    def get_all_taxa_names(self):
        """Wrapper method for the Root get_all_taxa_names method"""
        return self.files[0].get_all_taxa_names()

    def get_repeatpeps(self):
        """Wrapper method for the Root get_repeatpeps method"""
        return self.files[0].get_repeatpeps()

    # Leaf Wrapper methods ---------------------------------------------------------------------------------------
    def get_families_for_taxon(self, tax_id, curated_only=False, uncurated_only=False):
        """Returns accessions for all families directly associated with tax_id.
        Delegates to the root file's Lookup/ByTaxon, which covers all families."""
        return self.files[0].get_families_for_taxon(tax_id, curated_only, uncurated_only)

    def get_family_by_accession(self, accession, component=None):
        """Returns the family record for 'accession'.

        If 'component' is specified, only that component's files are searched.
        Otherwise all loaded files are searched (consensus files first, then hmm).
        """
        if component is not None:
            for leaf in self.components.get(component, {}).values():
                fam = leaf.get_family_by_accession(accession)
                if fam:
                    return fam
            return None

        # Search consensus components first (full metadata), then hmm
        search_order = (
            list(self.components[COMPONENT_CC].values()) +
            list(self.components[COMPONENT_UC].values()) +
            list(self.components[COMPONENT_CH].values()) +
            list(self.components[COMPONENT_UH].values())
        )
        for leaf in search_order:
            fam = leaf.get_family_by_accession(accession)
            if fam:
                return fam
        return None

    def get_family_by_accession_merged(self, accession):
        """Returns a Family built from consensus metadata merged with pHMM model data.

        Loads the consensus record for full metadata+sequence, then overlays
        the model field from the corresponding hmm file if available.
        Falls back to get_family_by_accession() if no consensus file is loaded.
        """
        # Try consensus first
        fam = self.get_family_by_accession(accession, component=COMPONENT_CC) or \
              self.get_family_by_accession(accession, component=COMPONENT_UC)
        if fam is None:
            return self.get_family_by_accession(accession)

        # Overlay model from hmm file if available
        hmm_fam = self.get_family_by_accession(accession, component=COMPONENT_CH) or \
                  self.get_family_by_accession(accession, component=COMPONENT_UH)
        if hmm_fam and hmm_fam.model:
            fam.model = hmm_fam.model
            if hmm_fam.taxa_thresholds:
                fam.taxa_thresholds = hmm_fam.taxa_thresholds
            if hmm_fam.general_cutoff:
                fam.general_cutoff = hmm_fam.general_cutoff
        return fam

    def get_family_by_name(self, name):
        """Wrapper method to call the Leaf get_family_by_name"""
        for fobj in self.files.values():
            try:
                fam = fobj.get_family_by_name(name)
                if fam:
                    return fam
            except Exception:
                pass
        return None

    def finalize(self):
        """Wrapper method to call the Leaf finalize"""
        for fobj in self.files.values():
            fobj.finalize()

    def filter_stages(self, accession, stages):
        """Wrapper method to call the Leaf filter_stages"""
        for fobj in self.files.values():
            fam = fobj.get_family_by_accession(accession)
            if fam:
                return fobj.filter_stages(accession, stages)

    def update_description(self, new_desc):
        """Wrapper method to call the Leaf update_description"""
        for file in self.files:
            self.files[file].update_description(new_desc)

    def check_unique(self, family):
        for fobj in self.files.values():
            if not fobj.check_unique(family):
                return False
        return True

    # File Utils
    def close(self):
        """Closes this FamDB instance, making further use invalid."""
        for fobj in self.files.values():
            fobj.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    # This method is here because famdb_data_loaders.py imports dfamorm, which is not available to users
    @staticmethod
    def read_embl_families(filename, lookup, header_cb=None):
        """
        This method is here because famdb_data_loaders.py imports dfamorm, which is not available to users

        Iterates over Family objects from the .embl file 'filename'. The format
        should match the output format of to_embl(), but this is not thoroughly
        tested.

        'lookup' should be a dictionary of Species names (in the EMBL file) to
        taxonomy IDs.

        If specified, 'header_cb' will be invoked with the contents of the
        header text at the top of the file before the iteration is complete.

        TODO: This mechanism is a bit awkward and should perhaps be reworked.
        """

        def set_family_code(family, code, value):
            """
            Sets an attribute on 'family' based on the EMBL line starting with 'code'.
            For codes corresponding to list attributes, values are appended.
            """
            if code == "ID":
                match = re.match(r"(\S*)", value)
                acc = match.group(1)
                acc = acc.rstrip(";")
                family.accession = acc
            elif code == "NM":
                family.name = value
            elif code == "DE":
                family.description = value
            elif code == "CC":
                # TODO: Consider only recognizing these after seeing "RepeatMasker Annotations"

                matches = re.match(r"\s*Type:\s*(\S+)", value)
                if matches:
                    family.repeat_type = matches.group(1).strip()

                matches = re.match(r"\s*SubType:\s*(\S+)", value)
                if matches:
                    family.repeat_subtype = matches.group(1).strip()

                matches = re.search(r"Species:\s*(.+)", value)
                if matches:
                    for spec in matches.group(1).split(","):
                        name = spec.strip()
                        if name:
                            tax_id = lookup.get(name)
                            if tax_id is not None:
                                family.clades += [tax_id]
                            else:
                                name = name.replace("[", "")
                                name = name.replace("]", "")
                                tax_id = lookup.get(name.lower())
                                if tax_id is not None:
                                    family.clades += [tax_id]
                                else:
                                    LOGGER.warning(
                                        f"Could not find taxon for '{name}' upper or lower: line={value}, and ID={family.accession}"
                                    )
                matches = re.search(r"SearchStages:\s*(\S+)", value)
                if matches:
                    family.search_stages = matches.group(1).strip()

                matches = re.search(r"BufferStages:\s*(\S+)", value)
                if matches:
                    family.buffer_stages = matches.group(1).strip()

                matches = re.search(r"Refineable", value)
                if matches:
                    family.refineable = True

        header = ""
        family = None
        in_header = True
        in_metadata = False

        nodes = lookup.values()

        with open(filename) as file:
            for line in file:
                if family is None:
                    # ID indicates start of metadata
                    if line.startswith("ID"):
                        family = Family()
                        family.clades = []
                        in_header = False
                        in_metadata = True
                    elif in_header:
                        matches = re.match(r"(CC)?\s*(.*)", line)
                        if line.startswith("XX"):
                            in_header = False
                        elif matches:
                            header_line = matches.group(2).rstrip("*").strip()
                            header += header_line + "\n"
                        else:
                            header += line

                if family is not None:
                    if in_metadata:
                        # SQ line indicates start of sequence
                        if line.startswith("SQ"):
                            in_metadata = False
                            family.consensus = ""

                        # Continuing metadata
                        else:
                            split = line.rstrip("\n").split(None, maxsplit=1)
                            if len(split) > 1:
                                code = split[0].strip()
                                value = split[1].strip()
                                set_family_code(family, code, value)

                    # '//' line indicates end of the sequence area
                    elif line.startswith("//"):
                        family.length = len(family.consensus)
                        keep = False
                        for clade in family.clades:
                            if clade in nodes:
                                LOGGER.debug(
                                    f"Including {family.accession} in taxa {clade} from {filename}"
                                )
                                keep = True
                        if keep:
                            yield family
                        family = None

                    # Part of the sequence area
                    else:
                        family.consensus += re.sub(r"[^A-Za-z]", "", line)

        # if header_cb:
        #     header_cb(header)