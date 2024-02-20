

def intarna_to_bplist(bplist_string):
    """Convert base pair list from intarna string format to python list.

    Arguments:
    bplist_string -- string with base pairs as for example returned by
                     IntaRNA. For example:
                    '(134,56):(135,55):(136,54):(137,53):(138,52):(139,51)'


    Returns:
    List of base pairs. Each base pair is represented as a tuple of the
    pairing positions. For example:
    [(134, 56), (135, 55), (136, 54), (137, 53), (138, 52), (139, 51)]

    """
    return [
        (
            int(item.split(",")[0].strip("(")),
            int(item.split(",")[1].strip(")")),
        )
        for item in bplist_string.split(":")
    ]


def get_string_representations(seq1, seq2, bp_list, id1="Seq1", id2="Seq2"):
    """Get interaction represented as a single string with three lines.

    Arguments:
    seq1 -- full sequence of first RNA (string)
    seq2 -- full sequence of second RNA (string)
    bp_list -- list of interacting base pairs (list of tuples; one based indices)
    id1 -- name of first RNA (string)
    id2 -- name of second RNA (string)

    Output:
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
