
import os
import subprocess
import pandas as pd
from io import StringIO
from typing import List, Optional, Union

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


def get_RRI_string_representations(
    seq1: str,
    seq2: str,
    bp_list: List[Tuple[int, int]],
    id1: str = "Seq1",
    id2: str = "Seq2"
) -> str:
    """
    Generate a string representation of RNA interaction, as a 3-line string.

    Args:
        seq1 (str): Full sequence of the first RNA.
        seq2 (str): Full sequence of the second RNA.
        bp_list (List[Tuple[int, int]]): List of interacting base pairs (1-based).
        id1 (str): Name/identifier for the first RNA.
        id2 (str): Name/identifier for the second RNA.

    Returns:
        str: Multi-line string representation of the interaction.

    Notes:
        The three-line string representation depicts the interaction between
        two RNAs within the interaction site. 
        - **First and Third Lines:** Represent the sequences of the two pairing
            RNAs, aligned with gaps to ensure pairing positions are properly matched. 
        - **Direction Annotation:** Sequence directions are explicitly marked with 5' and 3'. 
        - **Subsequence Annotation:** The indices of the first and last 
            nucleotides within the interaction site are provided in parentheses
            after the respective sequence IDs. These indices are 1-based.
        - **Second Line:** Displays the base pair interactions with vertical
            pipes (`|`) at aligned positions. Spaces indicate interior loops
            or bulges within the interaction site.

        Examples:

            5'-UACGGC-3' ArcZ[50:55]
            ||||||
            3'-AUGUCG-5' CyaR[34:29]

            5'-GAUUUCCUGGUGUAACGAAUUUUUUAAGUGC-3' DsrA[10:40]
            ||||||||  |||||||||||||  ||||||
            3'-CUAAAGGGGAACAUUGCUUAAAGU-UUUACG-5' rpoS[104:75]
    """

    # introduce gaps such that pairing sequence positions are aligned
    # and introduce pipes to mark pairing positions
    gapped_seq1, gapped_seq2, bps_as_string = "", "", ""

    for i in range(len(bp_list) - 1):
        len_a_frag =  - bp_list[i][0] + bp_list[i + 1][0]
        len_b_frag =    bp_list[i][1] - bp_list[i + 1][1]
        fragment_length = max(len_a_frag, len_b_frag)

        gapped_seq1 += (
            seq1[bp_list[i][0] - 1 : bp_list[i + 1][0] - 1]
            + (fragment_length - len_a_frag) * "-"
        )
        gapped_seq2 += (
            seq2[bp_list[i][1] - 1 : bp_list[i + 1][1] - 1 : -1]
            + (fragment_length - len_b_frag) * "-"
        )
        bps_as_string += "|" + " " * (fragment_length - 1)

    gapped_seq1 += seq1[bp_list[-1][0] - 1]
    gapped_seq2 += seq2[bp_list[-1][1] - 1]
    bps_as_string += "|"

    # annotate sequences
    gapped_seq1 = f"5'-{gapped_seq1}-3' {id1}[{bp_list[0][0]}:{bp_list[-1][0]}]"
    gapped_seq2 = f"3'-{gapped_seq2}-5' {id2}[{bp_list[0][1]}:{bp_list[-1][1]}]"
    bps_as_string = f"   {bps_as_string}"

    # unify length of line and return
    max_length = max(len(gapped_seq1), len(gapped_seq2))
    return "\n".join([
        gapped_seq1.ljust(max_length),
        bps_as_string.ljust(max_length),
        gapped_seq2.ljust(max_length)
    ])


def run_intarna(
    seq1: str,
    seq2: str,
    id1: str = "target",
    id2: str = "query",
    temperature: float = 37.0,
    intarna_args: Optional[List[str]] = None,
    out_file: Optional[str] = None,
    intarna_executable: str = "IntaRNA",
    outMode: str = "C",
    outCsvCols: str = "id1,id2,start1,end1,start2,end2,seq1,seq2,bpList,E,Etotal,"
                      "ED1,ED2,Pu1,Pu2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,"
                      "E_endR,E_hybrid,E_norm,E_add,P_E,hybridDPfull"
) -> Union[pd.DataFrame, str]:
    """
    Execute IntaRNA with specified parameters.

    Args:
        seq1 (str): Sequence of the first RNA.
        seq2 (str): Sequence of the second RNA.
        id1 (str): Identifier for the first RNA.
        id2 (str): Identifier for the second RNA.
        temperature (float): Temperature for the interaction prediction in Celsius.
        intarna_args (Optional[List[str]]): Additional arguments for IntaRNA.
        out_file (Optional[str]): Path to save the output (if provided).
        intarna_executable (str): Path or name of the IntaRNA executable.
        outMode (str): Output mode for IntaRNA.
        outCsvCols (str): Columns for CSV output (if outMode is "C").

    Returns:
        pd.DataFrame or str: Parsed DataFrame if outMode is "C"; otherwise, raw output.
    
    Notes:
        Only tested with csv output format.
    """
    if intarna_args is None:
        intarna_args = []

    args = [
        intarna_executable,
        "-t", seq1,
        "-q", seq2,
        "--tId", id1,
        "--qId", id2,
        "--temperature", str(temperature),
        "--outMode", outMode,
        "--outCsvCols="+ outCsvCols,
    ] + intarna_args

    if out_file:
        args.extend(["--out", out_file])

    result = subprocess.run(
        args,
        capture_output=True,
        text=True,
        universal_newlines=True,  # TODO: needed? pd needs stream anyway
        # stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE,
    )

    if result.returncode != 0:
        raise RuntimeError(f"IntaRNA failed with return code {result.returncode}: {result.stderr}")

    if outMode == "C":
        if out_file:
            return pd.read_csv(out_file, sep=";", comment="#")
        else:
            return pd.read_csv(StringIO(result.stdout), sep=";", comment="#")
    else:
        if out_file:
            with open(out_file, "r") as file:
                return file.read()
        return result.stdout