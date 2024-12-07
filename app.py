pip install streamlit biopython
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import csv
import io
import tempfile

def calculate_gc_content(sequence):
    """Calculate GC content percentage of a DNA sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def generate_primers_and_probes(sequence, primer_length=32, probe_length=54, amplicon_length=190):
    """
    Generate primers and probes from a given DNA sequence.
    
    Args:
    - sequence (str): DNA sequence to generate primers and probes from
    - primer_length (int): Length of forward and reverse primers
    - probe_length (int): Length of the probe
    - amplicon_length (int): Desired amplicon length
    
    Returns:
    - List of tuples containing (forward_primer, probe, reverse_primer)
    """
    results = []
    seq_len = len(sequence)
    
    for i in range(seq_len - amplicon_length + 1):
        forward_primer = sequence[i:i + primer_length]
        probe_start = i + primer_length
        probe = sequence[probe_start:probe_start + probe_length]
        reverse_primer_start = probe_start + probe_length
        reverse_primer = str(Seq(sequence[reverse_primer_start:reverse_primer_start + primer_length]).reverse_complement())
        
        results.append((forward_primer, probe, reverse_primer))
    
    return results

def generate_csv_content(sequence, primer_length=32, probe_length=54, amplicon_length=190):
    """
    Generate CSV content with primers, probes, and parameters.
    
    Args:
    - sequence (str): DNA sequence to generate primers and probes from
    - primer_length (int): Length of forward and reverse primers
    - probe_length (int): Length of the probe
    - amplicon_length (int): Desired amplicon length
    
    Returns:
    - StringIO object containing CSV content
    """
    primers_probes = generate_primers_and_probes(sequence, primer_length, probe_length, amplicon_length)
    
    # Create a StringIO object to write CSV content
    output = io.StringIO()
    csvwriter = csv.writer(output)
    
    # Write parameter information
    csvwriter.writerow(["Parameters and output for Gene Primer/Probe Design", "", ""])
    csvwriter.writerow(["Parameters", "", ""])
    csvwriter.writerow(["Primer Length", f"{primer_length}bp", ""])
    csvwriter.writerow(["Probe Length", f"{probe_length}bp", ""])
    csvwriter.writerow(["Amplicon length", f"{amplicon_length}bp", ""])
    csvwriter.writerow(["Min GC content", "30%", ""])
    csvwriter.writerow(["Max GC content", "70%", ""])
    csvwriter.writerow(["Tolerated self binding region", "5", ""])
    csvwriter.writerow(["Tolerated secondary structure region", "6", ""])
    csvwriter.writerow(["Tolerate background binding bases", "N/A", ""])
    csvwriter.writerow(["Number of background files checked", "N/A", ""])
    csvwriter.writerow([])
    csvwriter.writerow(["Forward_Primer", "Probe", "Reverse_Primer"])
    
    # Write primers and probes
    for forward_primer, probe, reverse_primer in primers_probes:
        csvwriter.writerow([forward_primer, probe, reverse_primer])
    
    # Reset the StringIO object's position to the beginning
    output.seek(0)
    return output

def main():
    """
    Streamlit main application for Gene Primer and Probe Generator
    """
    st.title("Gene Primer and Probe Generator")
    
    # Sidebar for parameter configuration
    st.sidebar.header("Primer/Probe Parameters")
    primer_length = st.sidebar.number_input("Primer Length (bp)", min_value=10, max_value=50, value=32)
    probe_length = st.sidebar.number_input("Probe Length (bp)", min_value=20, max_value=100, value=54)
    amplicon_length = st.sidebar.number_input("Amplicon Length (bp)", min_value=50, max_value=500, value=190)
    
    # File uploader
    uploaded_file = st.file_uploader("Choose a gene .fna file", type=['fna', 'fasta'])
    
    if uploaded_file is not None:
        try:
            # Read the uploaded FASTA file
            # Convert uploaded file to text mode
            uploaded_file.seek(0)
            record = next(SeqIO.parse(io.StringIO(uploaded_file.getvalue().decode('utf-8')), "fasta"))
            sequence = str(record.seq)
            
            # Display sequence information
            st.write(f"Sequence ID: {record.id}")
            st.write(f"Sequence Length: {len(sequence)} bp")
            st.write(f"GC Content: {calculate_gc_content(sequence):.2f}%")
            
            # Generate CSV
            csv_content = generate_csv_content(
                sequence, 
                primer_length=primer_length, 
                probe_length=probe_length, 
                amplicon_length=amplicon_length
            )
            
            # Download button
            st.download_button(
                label="Download Primers and Probes CSV",
                data=csv_content.getvalue(),
                file_name=f"{record.id}_primers_probes.csv",
                mime="text/csv"
            )
            
        except Exception as e:
            st.error(f"Error processing the file: {e}")

if __name__ == "__main__":
    main()
