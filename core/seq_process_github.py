import json
import os
import subprocess
import logging
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SequenceProcessingError(Exception):
    """Custom exception for sequence processing errors."""
    pass


def load_clustalw_path_from_config(config_path: str = "config.json") -> str:
    """
    Load the ClustalW executable path from a JSON config file.
    Raises SequenceProcessingError if the config or path is invalid.
    """
    if not os.path.exists(config_path):
        raise SequenceProcessingError(f"Config file '{config_path}' not found.")

    with open(config_path, "r") as f:
        config = json.load(f)

    clustalw_path = config.get("clustalw_path")
    if not clustalw_path or not os.path.isfile(clustalw_path):
        raise SequenceProcessingError("ClustalW executable path is missing or invalid in config.json.")

    return clustalw_path


def create_fasta_from_sequences(sequences, output_filepath: str, align: bool = True) -> str:
    """
    Create a FASTA file from a list of SequenceData objects, with optional ClustalW alignment.

    Parameters:
        sequences: List of SequenceData objects.
        output_filepath: Final FASTA output file path.
        align: Whether to align sequences using ClustalW.

    Returns:
        Path to the created FASTA file (aligned if `align` is True).
    """
    if not sequences:
        raise SequenceProcessingError("No sequences provided to create FASTA file.")

    fasta_records = []
    for seq_data in sequences:
        fasta_id = f"{seq_data.organism_name.replace(' ', '_')}_{seq_data.marker_code}_{seq_data.processid}"
        clean_sequence = "".join(filter(str.isalpha, seq_data.sequence.upper()))
        if not clean_sequence:
            logger.warning(f"Skipping empty or invalid sequence for ID {seq_data.processid}")
            continue

        record = SeqRecord(
            Seq(clean_sequence),
            id=fasta_id,
            name=seq_data.organism_name,
            description=f"{seq_data.organism_name} | {seq_data.marker_code} | BOLD ProcessID:{seq_data.processid}"
        )
        fasta_records.append(record)

    if not fasta_records:
        raise SequenceProcessingError("All sequences were invalid or empty.")

    output_dir = os.path.dirname(output_filepath)
    os.makedirs(output_dir, exist_ok=True)

    raw_fasta_path = output_filepath.replace(".fasta", "_unaligned.fasta")
    with open(raw_fasta_path, "w") as raw_handle:
        SeqIO.write(fasta_records, raw_handle, "fasta")
    logger.info(f"Raw FASTA file created: {raw_fasta_path}")

    if not align:
        return raw_fasta_path

    clustalw_path = load_clustalw_path_from_config()
    input_filename = os.path.basename(raw_fasta_path)
    clustalw_cmd = [clustalw_path, f"-infile={input_filename}"]

    logger.info(f"Running ClustalW: {' '.join(clustalw_cmd)} in directory '{output_dir}'")

    result = subprocess.run(
        clustalw_cmd,
        capture_output=True,
        text=True,
        cwd=output_dir
    )

    logger.debug("ClustalW STDOUT:\n" + result.stdout)
    logger.debug("ClustalW STDERR:\n" + result.stderr)

    if result.returncode != 0:
        raise SequenceProcessingError(f"ClustalW failed with return code {result.returncode}.\n{result.stderr}")

    aln_filename = os.path.splitext(input_filename)[0] + ".aln"
    aln_path = os.path.join(output_dir, aln_filename)

    if not os.path.isfile(aln_path):
        files_in_dir = "\n".join(os.listdir(output_dir))
        raise SequenceProcessingError(
            f"Alignment failed: .aln output file not found.\nFiles in '{output_dir}':\n{files_in_dir}"
        )

    alignment = AlignIO.read(aln_path, "clustal")
    with open(output_filepath, "w") as aligned_handle:
        AlignIO.write(alignment, aligned_handle, "fasta")

    logger.info(f"Aligned FASTA file created: {output_filepath}")
    return output_filepath
