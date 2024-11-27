

   def intarna_to_bplist(bplist_string: str, zero_based: bool = False) -> List[Tuple[int, int]]:
    """
    Convert a base pair list from IntaRNA string format to a Python list.

    Args:
        bplist_string (str): String with base pairs as returned by IntaRNA.
            Example: '(134,56):(135,55):(136,54):(137,53):(138,52):(139,51)'
        zero_based (bool): Convert 1-based nucleotide indices from IntaRNA to 0-based indices if True.

    Returns:
        List[Tuple[int, int]]: List of base pairs as tuples. Example:
        [(134, 56), (135, 55), ..., (139, 51)] or zero-based if `zero_based` is True.
    """

    offset = 1 if zero_based else 0
    return [
        (
            int(item.split(",")[0].strip("(")) - offset,
            int(item.split(",")[1].strip(")")) - offset,
        )
        for item in bplist_string.split(":")
    ]


def get_string_representations(seq1, seq2, bp_list, id1="Seq1", id2="Seq2"):
    """Get interaction represented as a single string with three lines.

    Args:
        seq1 (str): full sequence of first RNA
        seq2 (str): full sequence of second RNA
        bp_list (list): list of interacting base pairs in tuple. Note that indices should be one based
        id1 (str): name of first RNA
        id2 (str): name of second RNA

    Returns:
        Within the three line string representation, the first and third line
        repesenting the seqeunce of the two pairing RNAs within the
        interaction site. The sequences contain gaps such that the paring
        positions are aligned. The sequence directions are annotated with 5'
        and 3'. The subsequence is annotated after the sequence id by the
        (one based) index of the first and last nucleotide within the
        interaction site. Base pairs are marked by pipes in the
        corresponding positions within the second line. Interior loops and
        buldges within the interaction site correspond to spaces within the
        second line.

        Examples (missing tailing spaces):

        5'-UACGGC-3' ArcZ[50:55]
           ||||||
        3'-AUGUCG-5' CyaR[34:29]

        5'-GAUUUCCUGGUGUAACGAAUUUUUUAAGUGC-3' DsrA[10:40]
           ||||||||  |||||||||||||  ||||||
        3'-CUAAAGGGGAACAUUGCUUAAAGU-UUUACG-5' rpoS[104:75]
    """

    # introduce gaps such that pairing sequence positions are aligned
    # and introduce pipes to mark pairing positions

    gapped_seq1 = ""  # firs line
    gapped_seq2 = ""  # third line
    bps_as_string = ""  # second line
    for i in range(len(bp_list) - 1):
        len_a_frag = -bp_list[i][0] + bp_list[i + 1][0]
        len_b_frag = bp_list[i][1] - bp_list[i + 1][1]
        fragment_length = max(len_a_frag, len_b_frag)
        gapped_seq1 += (
            seq1[bp_list[i][0] - 1 : bp_list[i + 1][0] - 1]
            + (fragment_length - len_a_frag) * "-"
        )
        gapped_seq2 += (
            seq2[bp_list[i][1] - 1 : bp_list[i + 1][1] - 1 : -1]
            + (fragment_length - len_b_frag) * "-"
        )
        bps_as_string += "|" + (fragment_length - 1) * " "
    gapped_seq1 += seq1[bp_list[-1][0] - 1]
    gapped_seq2 += seq2[bp_list[-1][1] - 1]
    bps_as_string += "|"

    # annotate sequences
    gapped_seq1 = f"5'-{gapped_seq1}-3' {id1}[{bp_list[0][0]},{bp_list[-1][0]}]"
    gapped_seq2 = f"3'-{gapped_seq2}-5' {id2}[{bp_list[0][1]},{bp_list[-1][1]}]"
    bps_as_string = f"   {bps_as_string}    "

    # unify length of lines
    length = max([len(gapped_seq1), len(gapped_seq2)])
    gapped_seq1 = gapped_seq1.ljust(length)
    gapped_seq2 = gapped_seq2.ljust(length)
    bps_as_string = bps_as_string.ljust(length)

    # lines = [gapped_seq1, gapped_seq2, bps_as_string]
    # length = max([len(l) in lines])
    # lines = [l.ljust(length) for l in lines]
    # gapped_seq1, gapped_seq2, bps_as_string = lines

    return "\n".join([gapped_seq1, bps_as_string, gapped_seq2])

def run_intarna(
    seq1,
    seq2,
    id1="target",
    id2="query",
    temperature=37.0,
    intarna_args=[],
    out_file=None,
    intarna_executable="IntaRNA",
    outMode="C",
    outCsvCols="id1,id2,start1,end1,start2,end2,seq1,seq2,"
    + "bpList,E,Etotal,ED1,ED2,Pu1,Pu2,E_init,E_loops,E_dangleL,"
    + "E_dangleR,E_endL,E_endR,E_hybrid,E_norm,E_add,P_E,hybridDPfull",
):
    """Run IntaRNA.
        Only tested with csv output format.
        Provide following intarna arguments through function arguments
        and not through intarna_args variable:
        - --out  as out_file
        - --outMode as outMode
        - -t as seq1
        - -q as seq2
        - --tId  as id1
        - --qId as id2

    Args:
        seq1: sequence 1
        seq2: sequence 2
        id1: sequence 1 identifier
        id2: sequence 2 identifier
        temperature: temperature
        intarna_args: additional IntaRNA arguments
        out_file: path to store output
        intarna_executable: callable IntaRNA tool
        outMode: IntaRNA output mode
        outCsVCols: output csv format
    """
    # set up intarna arguments
    intarna_args = [
        intarna_executable,
        "-t",
        seq1,
        "-q",
        seq2,
        "--tId",
        id1,
        "--qId",
        id2,
        "--temperature",
        str(temperature),
        "--outMode",
        outMode,
        "--outCsvCols=" + outCsvCols,
    ] + intarna_args
    if out_file is not None:
        intarna_args.append("--out")
        intarna_args.append(out_file)

    # call intarna
    cp = subprocess.run(
        intarna_args,
        universal_newlines=True,  # TODO: needed? pd needs stream anyway
        # stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE,
        capture_output=True,
    )

    # prepare return format
    if cp.returncode != 0:
        print("IntaRNA returncode is " + str(cp.returncode))
        print(cp)

    # output
    if outMode == "C":
        if out_file is None:
            df = pd.read_csv(StringIO(cp.stdout), sep=";", comment="#")
        else:
            df = pd.read_csv(out_file, sep=";", comment="#")
        return df
    else:
        if out_file is None:
            return cp.stdout
        else:
            with open(out_file, "r") as out_handle:
                intarna_output = out_handle.read()
            return intarna_output