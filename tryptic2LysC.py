from __future__ import print_function
import os
import sys
import re
from collections import deque
import itertools
import argparse
import numpy as np
import pandas as pd

class Fasta(object):
    """
    not instantiated with file, in order to add multiple files, by setting and parsing iteratively.
    search for an AccessionNumber, Organism, DataBase, sequence, or partial sequence
    get proteotypic peptides for Organism
    write fasta of subset (e.g. AccessionNumbers, Organism, etc.)
    """
    def __init__(self):
        """
        an2aaseq_dict: key=AccessionNumber, val=AAseq
        org2an_dict: key=Organism, val=ListOfAN
        db2an_dict: key=DataBase, val=ListOfAN
        an2header_dict: key=AccessionNumber, val=header # original header
        :return: None
        """
        organism_regexes = [
            r"^>(dbj|emb|na|sp|tr|gb|ref)\|(.*?)\|.*?OS=(.*?)\s(?:(?:\w+=\w+\s?)+)(?:NCBI_TaxID=(\d+))?",
            r"^>([s|t][p|r])\|([\w]+).*?OS=(.*?) [A-Z][A-Z]=",
            r"^>gi\|\d+\|(\w+)\|(\w+\.\d).*\[(.*)\]",
            r"^>([^ ]+) ([^\[]+)"]
        self.my_regex = re.compile('|'.join(organism_regexes))
        self.an2aaseq_dict = {}
        self.org2an_dict = {}
        self.an2org_dict = {}
        self.db2an_dict = {}
        self.an2header_dict = {}

    def set_file(self, fn):
        self.fasta_file = fn

    def get_file(self):
        return self.fasta_file

    def parse_fasta(self, unique=True):
        """
        parse fasta-file by entering infos into dicts
        print summary information
        number of entries total, non-redundant aaseqs, entries per DataBase, entries per Organism
        :return: Str
        """
        for entry in self.yield_entry():
            header, aaseq = entry
            match = re.search(self.my_regex, header)
            if match: # at least one regex found something
                match_groups = filter(None, match.groups())
                if len(match_groups) == 2:
                    an, organism = match_groups
                    database = "na"
                elif len(match_groups) == 3:
                    database, an, organism = match_groups
                    database = database.strip()
                else:
                    print("Nothing could be extracted from the following fasta-entry: ", "\n", header, aaseq)
                    raise StopIteration
                an = an.strip()
                organism = organism.strip() # #!!! stop
            else:
                print("Nothing could be extracted from the following fasta-entry: ", "\n", header, aaseq)
                raise StopIteration
            # now fill dicts
            if not self.an2aaseq_dict.has_key(an):
                self.an2aaseq_dict[an] = aaseq
            else: # AccessionNumbers should be unique identifiers
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2aaseq_dict[an], "\n", header)
                    raise StopIteration
                else:
                    pass
                    # print("AN double entry: {}".format(an))
            if not self.org2an_dict.has_key(organism):
                self.org2an_dict[organism] = [an]
            else:
                self.org2an_dict[organism].append(an)

            if not self.an2org_dict.has_key(an):
                self.an2org_dict[an] = organism
            else:
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2org_dict[an], "\n", aaseq)
                    raise StopIteration
                else:
                    pass
            if not self.db2an_dict.has_key(database):
                self.db2an_dict[database] = [an]
            else:
                self.db2an_dict[database].append(an)

            if not self.an2header_dict.has_key(an):
                self.an2header_dict[an] = header
            else:
                if unique:
                    print("AccessionNumbers should be unique identifiers!")
                    print(an, "\n", self.an2aaseq_dict[an], "\n", aaseq)
                    raise StopIteration
                else:
                    pass

    def yield_entry(self):
        """
        generator that yields one entry of a fasta-file at a time,
        as tuple (header, aaseq)
        :return: tuple(Str, Str)
        """
        fh = open(self.get_file(), "r")
        aaseq = ""
        header = ""
        did_first = False
        for line in fh:
            if line[0] == ">":
                if did_first:
                    if len(aaseq) > 0:
                        yield(header, aaseq)
                else:
                    did_first = True
                header = line.rstrip()
                aaseq = ""
            else:
                aaseq += line.strip().upper()
        if len(aaseq) > 0:
            yield(header, aaseq)
        fh.close()

    def get_header_from_an(self, an):
        """
        return header of given AccessionNumber as String
        :param an: String
        """
        try:
            return self.an2header_dict[an]
        except:
            raise KeyError

    def get_aaseq_from_an(self, an):
        """
        return AAseq (protein sequence) of given AccessionNumber as String
        :param an: String
        :return: String
        """
        return self.an2aaseq_dict[an]

##### adapted from Pyteomics
##### https://bitbucket.org/levitsky/pyteomics
std_amino_acids = ['Q', 'W', 'E', 'R', 'T', 'Y', 'I', 'P', 'A', 'S',
                   'D', 'F', 'G', 'H', 'K', 'L', 'C', 'V', 'N', 'M']

expasy_rules = {
    'arg-c': r'R',
    'asp-n': r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1': r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2': r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3': r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4': r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5': r'(?<=[LW]EH)D',
    'caspase 6': r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7': r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8': r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9': r'(?<=LEH)D',
    'caspase 10': r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain': r'R',
    'cnbr': r'M',
    'enterokinase': r'(?<=[DE]{3})K',
    'factor xa': r'(?<=[AFGILTVM][DE]G)R',
    'formic acid': r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b': r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc': r'K',
    'ntcb': r'\w(?=C)',
    'pepsin ph1.3': r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                    r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin ph2.0': r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                    r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k': r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin': r'[^DE](?=[AFILMV])',
    'thrombin': r'((?<=G)R(?=G))|'
                r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin': r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
}
"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the
`PeptideCutter tool
<http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
at Expasy.
"""

"""modX labels for the 20 standard amino acids."""

std_nterm = 'H-'
"""modX label for the unmodified N-terminus."""

std_cterm = '-OH'
"""modX label for the unmodified C-terminus."""

std_labels = std_amino_acids + [std_nterm, std_cterm]
"""modX labels for the standard amino acids and unmodified termini."""

def is_term_mod(label):
    """Check if `label` corresponds to a terminal modification.

    Parameters
    ----------
    label : str

    Returns
    -------
    out : bool
    """
    return label[0] == '-' or label[-1] == '-'

def length(sequence, **kwargs):
    """Calculate the number of amino acid residues in a polypeptide
    written in modX notation.

    Parameters
    ----------
    sequence : str or list or dict
        A string with a polypeptide sequence, a list with a parsed sequence or
        a dict of amino acid composition.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Examples
    --------
    >>> length('PEPTIDE')
    7
    >>> length('H-PEPTIDE-OH')
    7
    """
    if not sequence: return 0

    if isinstance(sequence, str) or isinstance(sequence, list):
        if isinstance(sequence, str):
            parsed_sequence = parse(sequence, **kwargs)
        else:
            parsed_sequence = sequence
        num_term_groups = 0
        if is_term_mod(parsed_sequence[0]):
            num_term_groups += 1
        if is_term_mod(parsed_sequence[-1]):
            num_term_groups += 1
        return len(parsed_sequence) - num_term_groups
    elif isinstance(sequence, dict):
        return sum(amount for aa, amount in sequence.items()
                   if not is_term_mod(aa))

    print('Unsupported type of sequence.')


_modX_sequence = re.compile(r'^([^-]+-)?((?:[a-z]*[A-Z])+)(-[^-]+)?$')
_modX_group = re.compile(r'[a-z]*[A-Z]')
_modX_split = re.compile(r'([a-z]*)([A-Z])')

def parse(sequence, show_unmodified_termini=False, split=False,
          allow_unknown_modifications=False,
          **kwargs):
    """Parse a sequence string written in modX notation into a list of
    labels or (if `split` argument is :py:const:`True`) into a list of
    tuples representing amino acid residues and their modifications.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    show_unmodified_termini : bool, optional
        If :py:const:`True` then the unmodified N- and C-termini are explicitly
        shown in the returned list. Default value is :py:const:`False`.
    split : bool, optional
        If :py:const:`True` then the result will be a list of tuples with 1 to 4
        elements: terminal modification, modification, residue. Default value is
        :py:const:`False`.
    allow_unknown_modifications : bool, optional
        If :py:const:`True` then do not raise an exception when an unknown
        modification of a known amino acid residue is found in the sequence.
        This also includes terminal groups.
        Default value is :py:const:`False`.

        .. note::
            Since version 2.5, this parameter has effect only if `labels`
            are provided.
    labels : container, optional
        A container of allowed labels for amino acids,
        modifications and terminal modifications.
        If not provided, no checks will be done.
        Separate labels for modifications (such as 'p' or 'ox')
        can be supplied, which means they are applicable to all residues.

        .. warning::
            If `show_unmodified_termini` is set to :py:const:`True`, standard
            terminal groups need to be present in `labels`.

        .. warning::
            Avoid using sequences with only one terminal group, as they are
            ambiguous. If you provide one, `labels` (or :py:const:`std_labels`)
            will be used to resolve the ambiguity.

    Returns
    -------
    out : list
        List of tuples with labels of modifications and amino acid residues.

    Examples
    --------
    >>> parse('PEPTIDE', split=True)
    [('P',), ('E',), ('P',), ('T',), ('I',), ('D',), ('E',)]
    >>> parse('H-PEPTIDE')
    ['P', 'E', 'P', 'T', 'I', 'D', 'E']
    >>> parse('PEPTIDE', show_unmodified_termini=True)
    ['H-', 'P', 'E', 'P', 'T', 'I', 'D', 'E', '-OH']
    >>> parse('TEpSToxM', labels=std_labels + ['pS', 'oxM'])
    ['T', 'E', 'pS', 'T', 'oxM']
    >>> parse('zPEPzTIDzE', True, True, labels=std_labels+['z'])
    [('H-', 'z', 'P'), ('E',), ('P',), ('z', 'T'), ('I',), ('D',), ('z', 'E', '-OH')]
    """
    sequence = str(sequence)

    try:
        n, body, c = re.match(_modX_sequence, sequence).groups()
    except AttributeError:
        print('Not a valid modX sequence: ' + sequence)

    # Check for allowed labels, if they were explicitly given
    labels = kwargs.get('labels')
    # labels help save the day when only one terminal group is given
    if c is None and n is not None:
        if labels is None:
            labels = std_labels
        # we can try to resolve the ambiguity
        if n != std_nterm and n not in labels:
            # n is the body then
            c = '-' + body
            body = n[:-1]
            n = None

    # Actual parsing
    if split:
        parsed_sequence = [g if g[0] else (g[1],) for g in re.findall(
            _modX_split, body)]
    else:
        parsed_sequence = re.findall(_modX_group, body)
    nterm, cterm = (n or std_nterm), (c or std_cterm)

    # Check against `labels` if given
    if labels is not None:
        labels = set(labels)
        for term, std_term in zip([n, c], [std_nterm, std_cterm]):
            if term and term not in labels and not allow_unknown_modifications:
                print(
                    'Unknown label: {}'.format(term))
        for group in parsed_sequence:
            if split:
                mod, X = group if len(group) == 2 else ('', group[0])
            else:
                mod, X = re.match(_modX_split, group).groups()
            if ((not mod) and X not in labels) or not ((mod + X in labels) or (
                            X in labels and (
                                    mod in labels or allow_unknown_modifications))):
                print(
                    'Unknown label: {}'.format(group))

    # Append terminal labels
    if show_unmodified_termini or nterm != std_nterm:
        if split:
            parsed_sequence[0] = (nterm,) + parsed_sequence[0]
        else:
            parsed_sequence.insert(0, nterm)
    if show_unmodified_termini or cterm != std_cterm:
        if split:
            parsed_sequence[-1] = parsed_sequence[-1] + (cterm,)
        else:
            parsed_sequence.append(cterm)

    return parsed_sequence

def cleave(sequence, rule, missed_cleavages=0, min_length=None, **kwargs):
    """Cleaves a polypeptide sequence using a given rule.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str or compiled regex
        A regular expression describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose C-terminal
        bond is to be cleaved. All additional requirements should be specified
        using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
    missed_cleavages : int, optional
        Maximum number of allowed missed cleavages. Defaults to 0.
    min_length : int or None, optional
        Minimum peptide length. Defaults to :py:const:`None`.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Returns
    -------
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0, labels='ABK') == {'AK', 'BK'}
    True
    >>> cleave('GKGKYKCK', expasy_rules['trypsin'], 2) == \
    {'CK', 'GKYK', 'YKCK', 'GKGK', 'GKYKCK', 'GK', 'GKGKYK', 'YK'}
    True

    """
    sequence = sequence.replace("*", "")
    return sorted(set(_cleave(sequence, rule, missed_cleavages, min_length, **kwargs)))

def _cleave(sequence, rule, missed_cleavages=2, min_length=7, **kwargs):
    """Like :py:func:`cleave`, but the result is a list. Refer to
    :py:func:`cleave` for explanation of parameters.
    """
    peptides = []
    cleavage_sites = deque([0], maxlen=missed_cleavages + 2)
    for i in itertools.chain(map(lambda x: x.end(), re.finditer(rule, sequence)), [None]):
        cleavage_sites.append(i)
        for j in range(len(cleavage_sites) - 1):
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or length(seq, **kwargs) >= min_length:
                    peptides.append(seq)
    return peptides

def clean_an(an):
    return an.replace("CON__", "").split("-")[0] #.replace("REV__", "")

def num_missed_cleavages_lysc(peptide_sequence):
    try:
        mc = peptide_sequence.count("K") - 1
    except AttributeError:
        return np.nan
    if mc < 0:
        return 0
    else:
        return mc

def get_start_stop_position_of_sequence(df, fa):
    an = df["Leading razor protein"]
    try:
        aaseq = fa.get_aaseq_from_an(an)
    except KeyError:
        try:
            aaseq = fa.get_aaseq_from_an(clean_an(an))
        except KeyError:
            return np.nan, np.nan
    peptide = df["Sequence"]
    try:
        start = aaseq.index(peptide)
    except ValueError:
        return np.nan, np.nan
    stop = start + len(peptide)
    return start, stop

def get_least_mc_hit(start_end_mc_list):
    start_end_mc_list.sort(key=lambda x: x[-1], reverse=False)
    return start_end_mc_list[0]

def get_start_stop_missed_cleavage_LysC_v2(df, AN2LysC_dict):
    start_tryp = df["Start_Tryp"]
    end_tryp = df["End_Tryp"]
    if np.isnan(start_tryp) or np.isnan(end_tryp):
        return np.nan, np.nan, np.nan
    an = df["Leading razor protein"]
    try:
        start_end_mc_list = AN2LysC_dict[an]
    except KeyError:
        try:
            start_end_mc_list = AN2LysC_dict[clean_an(an)]
        except KeyError:
            return np.nan, np.nan, np.nan
    start_end_mc_candidates = []
    for start_end_mc in start_end_mc_list:
        start_lysc, end_lysc, mc_lysc = start_end_mc
        if start_tryp >= start_lysc and end_tryp <= end_lysc:
            start_end_mc_candidates.append(start_end_mc)
    if not start_end_mc_candidates:
        print(an)
        return np.nan, np.nan, np.nan
    else: # if any hits, get the one with least mc
        return get_least_mc_hit(start_end_mc_candidates)

def get_AN2LysC_dict_v2(fa, ans_list):
    #### LysC precalculation
    # AN2LysC_dict = key: AN, val: ListOfTuple(Start, End position, missed_cleavages)
    AN2LysC_dict = {}
    for an in ans_list:
        try:
            aaseq = fa.get_aaseq_from_an(an)
        except KeyError:
            try:
                aaseq = fa.get_aaseq_from_an(clean_an(an))
            except KeyError:
                continue
        # peptide_list = lysc_digestion(aaseq, miss_cleavage=6, min_length=7)
        peptide_list = cleave(aaseq, expasy_rules["lysc"], missed_cleavages=6, min_length=7)
        start_end_mc_list = []
        for peptide in peptide_list:
            start = aaseq.index(peptide)
            end = start + len(peptide)
            missed_cleavages = num_missed_cleavages_lysc(peptide)
            start_end_mc = (start, end, missed_cleavages)
            start_end_mc_list.append(start_end_mc)
        if not an in AN2LysC_dict:
            AN2LysC_dict[an] = start_end_mc_list
        else:
            print(an)
    return AN2LysC_dict

def get_lysc_peptide(df, fa):
    start_lysc = df["Start_LysC"]
    end_lysc = df["End_LysC"]
    if np.isnan(start_lysc) or np.isnan(end_lysc):
        return np.nan
    an = df["Leading razor protein"]
    try:
        return fa.get_aaseq_from_an(an)[int(start_lysc):int(end_lysc)]
    except KeyError:
        try:
            return fa.get_aaseq_from_an(clean_an(an))[int(start_lysc):int(end_lysc)]
        except KeyError:
            return np.nan

def grep_gene_symbol(header):
    try:
        return header.split("|")[2].split(" ")[0]
    except AttributeError:
        return np.nan

def get_header_from_an(an, fa):
    try:
        return fa.get_header_from_an(an)
    except KeyError:
        return np.nan

def get_bufferedPeptide_from_aaseq_and_pos(aaseq, position, window_size):
    """
    e.g. VQPSSSLSARMMSGSRGSSLNDTYHSRDSSF or __________MSTGLRYKSKLATPEDKQDID or LETSEDDSEEDGDHNRTFDI___________
    param aaseq: String(Amino Acid sequence)
    param position: Integer(starting count from 0, NOT m)
    param window_size: positive Integer, bi-directional number of AAs to include
    """
    len_aaseq = len(aaseq) - 1
    if position > len_aaseq:
        raise IndexError
    start_index = position - window_size
    stop_index = position + window_size + 1
    if start_index >= 0 and stop_index <= len_aaseq:
        return aaseq[start_index : stop_index]
    elif start_index < 0:
        return "_" * abs(start_index) + aaseq[0 : stop_index]
    elif stop_index > len_aaseq:
        return aaseq[start_index : len_aaseq + 1] + "_" * (stop_index - len_aaseq - 1)

def LysCPepMod_ResMod_ResNum(df, aaseq_cut_length, mod):
    Start_LysC = int(df["Start_LysC"])
    End_LysC = int(df["End_LysC"]) - 1
    Residue_number_list = range(Start_LysC, End_LysC, aaseq_cut_length)
    # include the very last AA
    if not End_LysC in Residue_number_list:
        Residue_number_list.append(End_LysC)
    # Residue_number_list: positions within protein-aaseq with modification
    shift = abs(0 - Residue_number_list[0])
    Residue_number_list_shift = [ele - shift for ele in Residue_number_list]
    aaseq = df["LysC_pep"]
    aaseq_list = list(aaseq)
    Residue_mod_list = []
    for res in Residue_number_list_shift:
        Residue_mod_list.append(aaseq[res])
        aaseq_list[res] = aaseq_list[res] + mod
    aaseq_mod = "".join(aaseq_list)
    return aaseq_mod, Residue_mod_list, Residue_number_list

def melt_via_iterrows(df, fa, modification_name):
    modseq_list = []
    protmatch_list = []
    protname_list = []
    modname_list = []
    residue_list = []
    residuenumber_list = []
    preceedingAA_list = []
    trailingAA_list = []
    seqwindow_list = []
    protdescr_list = []
    database_list = []
    for index_, ser in df.iterrows():
        Residue_number_list = ser["Residue_number_list"]
        # LysC_pep = ser['LysC_pep']
        Start_LysC = int(ser['Start_LysC'])
        End_LysC = int(ser['End_LysC']) - 1
        protmatch = ser["PROTEIN MATCH / IPI ID"]
        an = protmatch
        protname = ser["PROTEIN NAME / GENE SYMBOL"]
        protein_description = ser["Protein Description"]
        database = ser["Database"]
        modseq = ser["LysC_pep_mod"]
        Residue_mod_list = ser["Residue_mod_list"]
        try:
            aaseq = fa.get_aaseq_from_an(an)
        except KeyError:
            continue
        if not Residue_number_list:
            continue
        if Start_LysC == 0:
            preceeding_aa = "-"
        else:
            preceeding_aa = aaseq[Start_LysC-1]
        if End_LysC == (len(aaseq) - 1):
            trailing_aa = "-"
        else:
            trailing_aa = aaseq[End_LysC+1]
        for index_, res_number in enumerate(Residue_number_list):
            # Residue_number_list and residue_list need to be of equal length
            # Residue_number_list is the index of AA in the protein
            # residue_list is the AA within LysC peptide that has been modified
            residue = Residue_mod_list[index_]
            modseq_list.append(modseq)
            protmatch_list.append(protmatch)
            protname_list.append(protname)
            modname_list.append(modification_name)
            residue_list.append(residue)
            residuenumber_list.append(res_number)
            preceedingAA_list.append(preceeding_aa)
            trailingAA_list.append(trailing_aa)
            if res_number == len(aaseq):
            #### "End_LysC" is for python slicing, but is NOT the index position (counting from 0) of the string
            #### therefore -1
                res_number -= 1
            buff_pep = get_bufferedPeptide_from_aaseq_and_pos(aaseq, res_number, window_size=6)
            seqwindow_list.append(buff_pep)
            protdescr_list.append(protein_description)
            database_list.append(database)

    dict2df = {'ModifiedSequence - Evidence Top Mascot Scoring Peptide': modseq_list,
               'PROTEIN MATCH / IPI ID': protmatch_list,
               'PROTEIN NAME / GENE SYMBOL': protname_list,
               'MODIFICATION NAME': modname_list,
               'RESIDUE': residue_list,
               'RESIDUE NUMBER': residuenumber_list,
               'PRECEEDING AA': preceedingAA_list,
               'TRAILING AA': trailingAA_list,
               'SEQUENCE WINDOW (Modification +/-6 AA)': seqwindow_list,
               'Protein Description': protdescr_list,
               'Database': database_list}
    return pd.DataFrame(dict2df)

def run(fn_fasta, fn_mq, fn_out_lysc, fn_out_ptm, aaseq_cut_length, modification_name):
    fa = Fasta()
    fa.set_file(fn_fasta)
    fa.parse_fasta()
    df = pd.read_csv(fn_mq, sep='\t')
    df.dropna(axis=0, how="all", inplace=True) # remove empty rows

    # AN2Tryp_dict = get_AN2Tryp_dict(fa)
    ans_list = df["Leading razor protein"].unique().tolist()
    AN2LysC_dict = get_AN2LysC_dict_v2(fa, ans_list) # changed method to pyteomics cleave
    df["Start_Tryp"], df["End_Tryp"] = zip(*df.apply(get_start_stop_position_of_sequence, args=(fa, ), axis=1))
    df["Start_LysC"], df["End_LysC"], df["missed_cleavages_LysC"] = zip(*df.apply(get_start_stop_missed_cleavage_LysC_v2, axis=1, args=(AN2LysC_dict,)))

    # sanity check: all peptides were successfully mapped (excluding reverse and contaminants)
    cond_rev = df["Reverse"].notnull()
    cond_cont = df["Potential contaminant"].notnull()
    cond_rev_cont = cond_rev | cond_cont
    assert sum(df.loc[-cond_rev_cont, "Start_Tryp"].isnull()) == 0
    assert sum(df.loc[-cond_rev_cont, "Start_LysC"].isnull()) == 0
    df["LysC_pep"] = df.apply(get_lysc_peptide, args=(fa, ), axis=1)
    df_out = df.copy()
    df_out["Start_Tryp"] += 1
    df_out["Start_LysC"] += 1
    df_out = df_out.fillna("")
    df_out.to_csv(fn_out_lysc, sep='\t', header=True, index=False)
    df = df[-cond_rev_cont]

####### PTM finder analogue file
    dfx = pd.DataFrame()
    dfx["LysC_pep"] = df.apply(get_lysc_peptide, args=(fa, ), axis=1)
    dfx["PROTEIN MATCH / IPI ID"] = df["Leading razor protein"]
    dfx["Start_LysC"] = df["Start_LysC"]
    dfx["End_LysC"] = df["End_LysC"]
    dfx["Database"] = os.path.basename(fn_fasta)
    dfx['Protein Description'] = df["Leading razor protein"].apply(get_header_from_an, args=(fa, ))
    dfx["PROTEIN NAME / GENE SYMBOL"] = dfx["Protein Description"].apply(grep_gene_symbol)
    dfx['MODIFICATION NAME'] = modification_name
    dfx["LysC_pep_mod"], dfx["Residue_mod_list"], dfx["Residue_number_list"] = zip(*dfx.apply(LysCPepMod_ResMod_ResNum, args=(aaseq_cut_length, modification_name, ), axis=1))
    dfm = melt_via_iterrows(dfx, fa, modification_name)
    COLS_PTMFINDER = ['ModifiedSequence - Evidence Top Mascot Scoring Peptide',
                      'PROTEIN MATCH / IPI ID', 'PROTEIN NAME / GENE SYMBOL',
                      'MODIFICATION NAME', 'RESIDUE', 'RESIDUE NUMBER',
                      'PRECEEDING AA', 'TRAILING AA',
                      'SEQUENCE WINDOW (Modification +/-6 AA)',
                      'Protein Description', 'Database']
    dfm = dfm[COLS_PTMFINDER]

    dfm["RESIDUE NUMBER"] += 1
    dfm.sort_values(["PROTEIN NAME / GENE SYMBOL", "PROTEIN MATCH / IPI ID", "RESIDUE NUMBER"],
                     ascending=[True, True, True], inplace=True)
    dfm.drop_duplicates(inplace=True)
    dfm.to_csv(fn_out_ptm, sep='\t', header=True, index=False)
    print("")
    print("\t Finished processing :)")
    print("#" * 80)

def error_(parser):
    sys.stderr.write("The arguments passed are invalid.\nPlease check the input parameters.\n\n")
    parser.print_help()
    sys.exit(2)

if __name__ == "__main__":
    cmd_line = True
    if cmd_line:
        parser = argparse.ArgumentParser()

        parser.add_argument("-fasta", "--fn_fasta", help="FASTA file absolute path", type=str)
                            # default=r"absolute/path/to/directory/filename.ending")

        parser.add_argument("-mq", "--fn_maxquant", help="MaxQuant peptides file absolute path", type=str)

        parser.add_argument("-lysc", "--fn_out_lysc", help="LysC output file absolute path", type=str,
                            default=None)

        parser.add_argument("-ptm", "--fn_out_ptm", help="PTM-finder-like output file absolute path", type=str,
                            default=None)

        parser.add_argument("-cut", "--aaseq_cut_length", help="Peptide length, number of amino acids to split LysC-peptides into (default=20)", type=int,
                            default=20)

        parser.add_argument("-mod", "--modification_name",
                            help="name of the amino acid modification (default='(ph)'",
                            type=str, default="(ph)")

        args = parser.parse_args()

        for arg in sorted(vars(args)):
            if getattr(args, arg) is None:
                error_(parser)

        print("#" * 80)
        for arg in sorted(vars(args)):
            print(arg, ": ", getattr(args, arg))

        fn_mq = args.fn_maxquant
        if not args.fn_out_lysc:
            fn_out_lysc = fn_mq.replace(".txt", "_LysC.txt")
        else:
            fn_out_lysc = args.fn_out_lysc
        if not args.fn_out_ptm:
            fn_out_ptm = fn_mq.replace(".txt", "_forPerseus.txt")
        else:
            fn_out_ptm = args.fn_out_ptm

        run(args.fn_fasta, fn_mq, fn_out_lysc, fn_out_ptm, args.aaseq_cut_length, args.modification_name)

    else:
        fn_fasta = r"absolute/path/to/directory/filename.ending"
        fn_mq = r"absolute/path/to/directory/filename.ending"
        fn_out_lysc = r"absolute/path/to/directory/filename.ending"
        fn_out_ptm = r"absolute/path/to/directory/filename.ending"
        aaseq_cut_length = 20
        modification_name = "(ph)"
        run(fn_fasta, fn_mq, fn_out_lysc, fn_out_ptm, aaseq_cut_length, modification_name)

