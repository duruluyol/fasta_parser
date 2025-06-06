import os
import random
import logging
import requests
from io import StringIO
from Bio import SeqIO

from core.sequence_processing import create_fasta_from_sequences, SequenceProcessingError

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

class BOLDAPIError(Exception):
    pass

class SequenceData:
    def __init__(self, processid, organism_name, marker_code, sequence, query_name):
        self.processid = processid
        self.organism_name = organism_name
        self.marker_code = marker_code
        self.sequence = sequence
        self.query_name = query_name

    def __repr__(self):
        return (f"SequenceData(processid='{self.processid}', "
                f"organism='{self.organism_name}', "
                f"marker='{self.marker_code}', "
                f"query_name='{self.query_name}', "
                f"sequence='{self.sequence[:30]}...')")

BOLD_API_URL = "http://www.boldsystems.org/index.php/API_Public"

def fetch_sequences_by_organism_and_marker(organism: str, marker_code: str, limit: int = 10) -> list[SequenceData]:
    params = {
        'taxon': organism,
        'marker': marker_code,
    }
    url = f"{BOLD_API_URL}/sequence"

    logging.debug(f"Requesting URL: {url} with params: {params}")

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise BOLDAPIError(f"Failed to connect to BOLD Systems API: {e}")

    if not response.text.strip():
        logging.warning("BOLD Systems API returned an empty or whitespace-only response.")
        return []

    if response.text.strip().lower().startswith("<!doctype html>") or \
       response.text.strip().lower().startswith("<html"):
        raise BOLDAPIError("BOLD Systems API returned an HTML page instead of FASTA data.")

    sequences = []
    try:
        for record in SeqIO.parse(StringIO(response.text), "fasta"):
            parts = record.id.split('|')
            processid = parts[0] if len(parts) > 0 else "N/A"
            organism_name = parts[1] if len(parts) > 1 else "N/A"
            marker_code_found = parts[2] if len(parts) > 2 else "N/A"

            sequences.append(
                SequenceData(
                    processid=processid,
                    organism_name=organism_name,
                    marker_code=marker_code_found,
                    sequence=str(record.seq),
                    query_name=organism
                )
            )
    except Exception as e:
        raise BOLDAPIError(f"Failed to parse FASTA response from BOLD Systems API. Error: {e}")

    if len(sequences) > limit:
        sequences = random.sample(sequences, limit)

    return sequences

if __name__ == "__main__":
    logging.info("--- BOLD Systems DNA Sequence Fetcher ---")

    organism_input = input("Enter organism names separated by commas (e.g. Homo sapiens, Mus musculus):\n> ")
    organisms_to_fetch = [org.strip() for org in organism_input.split(',') if org.strip()]
    common_barcode_type = input("Enter the barcode DNA type (e.g. COI-5P, 16S, ITS):\n> ").strip()

    try:
        sequences_per_organism_limit = int(input("Enter the number of sequences to fetch per organism:\n> ").strip())
    except ValueError:
        logging.warning("Invalid number entered. Using default value of 5.")
        sequences_per_organism_limit = 5

    all_fetched_sequences = []

    output_base_dir = os.path.join("data", "temp_fasta")
    combined_fasta_filename = f"combined_sequences_for_{common_barcode_type}.fasta"
    combined_fasta_path = os.path.join(output_base_dir, combined_fasta_filename)

    os.makedirs(output_base_dir, exist_ok=True)
    logging.info(f"Output directory ensured: {output_base_dir}")

    for organism in organisms_to_fetch:
        logging.info(f"Fetching {sequences_per_organism_limit} '{common_barcode_type}' sequences for: {organism}")
        try:
            current_organism_sequences = fetch_sequences_by_organism_and_marker(
                organism=organism,
                marker_code=common_barcode_type,
                limit=sequences_per_organism_limit
            )
            if current_organism_sequences:
                logging.info(f"  Successfully fetched {len(current_organism_sequences)} sequences for {organism}.")
                all_fetched_sequences.extend(current_organism_sequences)
            else:
                logging.warning(f"  No sequences found for {organism}.")
        except BOLDAPIError as e:
            logging.error(f"  Error fetching sequences for {organism}: {e}")
        except Exception as e:
            logging.exception(f"  An unexpected error occurred for {organism}.")

    if all_fetched_sequences:
        logging.info(f"--- Total sequences fetched: {len(all_fetched_sequences)} ---")
        logging.info(f"Creating combined FASTA file at: {combined_fasta_path}")
        try:
            final_fasta_path = create_fasta_from_sequences(all_fetched_sequences, combined_fasta_path)
            logging.info(f"FASTA file created successfully at: {final_fasta_path}")
        except SequenceProcessingError as e:
            logging.error(f"Error creating combined FASTA file: {e}")
        except Exception as e:
            logging.exception(f"Unexpected error during FASTA file creation.")
    else:
        logging.warning("No sequences were fetched. No FASTA file was created.")

    logging.info("--- Process Finished ---")
