from pymmseqs.commands import easy_cluster
import os

def run_easy_cluster(input_fasta, tmp_dir, output_dir, min_cluster_size=1):
    """
    Run MMseqs2 easy-cluster on a FASTA file with progress messages.
    
    Parameters:
        input_fasta (str): Path to input FASTA file.
        tmp_dir (str): Path to temporary directory.
        output_dir (str): Path to output directory.
        min_seq_id (float): Minimum sequence identity threshold.
        coverage (float): Coverage threshold.
    """
    
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("[1/2] Running MMseqs2 easy-cluster...")
    aa_cluster = easy_cluster(input_fasta, output_dir, tmp_dir, min_seq_id=0.7)

    print("[2/2] Saving cluster representatives to FASTA...")
    fasta_path = os.path.join(output_dir, "cluster_reps.fasta")
    with open(fasta_path, "w") as f:
        for cluster in aa_cluster.to_gen():
            if len(cluster["members"]) >= min_cluster_size:
                print(f"Representative sequence of a cluster: {cluster['rep']}")
                print(cluster.keys())
                rep_seq = cluster["rep"]
                rep_id = cluster["members"][0] # usually the sequence ID
                f.write(f">{rep_id['header']}\n{rep_id['sequence']}\n")

    print(f"Done! FASTA file written to: {fasta_path}")


if __name__ == "__main__":
    # Example usage
    run_easy_cluster(
        input_fasta="aa-uniref50.fasta",
        # input_fasta="mixed_split_meltome_flip.fasta",
        tmp_dir="/home/elias/mmseqs_tmp",
        output_dir="/home/elias/clusters",
        min_cluster_size=1
    )