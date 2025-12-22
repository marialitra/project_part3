import numpy as np
import sys

def read_embeddings(dat_file, num_proteins=10):
    """
    Read embeddings from binary .dat file with embedded protein IDs.
    
    Args:
        dat_file: Path to the .dat file (e.g., output/vectors_with_ids.dat)
        num_proteins: Number of proteins to display (default: 10)
    """
    try:
        # Configuration (must match protein_embed_danai.py)
        EMBED_DIM = 320
        MAX_ID_LEN = 50
        total_cols = EMBED_DIM + MAX_ID_LEN
        
        # Load memmap
        data = np.memmap(dat_file, dtype='float32', mode='r')
        num_rows = len(data) // total_cols
        data = data.reshape(num_rows, total_cols)
        
        # Extract embeddings and IDs
        embeddings = data[:, :EMBED_DIM]
        id_arrays = data[:, EMBED_DIM:]
        
        # Decode protein IDs
        protein_ids = []
        for id_arr in id_arrays:
            # Convert float32 back to bytes
            byte_arr = id_arr.astype(np.uint8).tobytes()
            # Decode and strip null bytes
            protein_id = byte_arr.decode('ascii', errors='replace').rstrip('\x00')
            protein_ids.append(protein_id)
        
        print(f"[INFO] Loaded {len(protein_ids)} proteins")
        print(f"[INFO] Embedding shape per protein: {embeddings.shape[1]}")
        print(f"[INFO] Total embedding array shape: {embeddings.shape}")
        
        # Display first x proteins
        print(f"\n[INFO] First {min(num_proteins, len(protein_ids))} proteins:")
        for i in range(min(num_proteins, len(protein_ids))):
            print(f"\nProtein: {protein_ids[i]}")
            print(f"Embedding (first 10 values): {embeddings[i][:10]}")
            print(f"Embedding stats - Min: {embeddings[i].min():.4f}, Max: {embeddings[i].max():.4f}, Mean: {embeddings[i].mean():.4f}")
        
    except FileNotFoundError as e:
        print(f"[ERROR] File not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Failed to read file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    # Default values
    dat_file = "output/vectors_with_ids.dat"
    num_proteins = 10
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        dat_file = sys.argv[1]
    if len(sys.argv) > 2:
        num_proteins = int(sys.argv[2])
    
    read_embeddings(dat_file, num_proteins)
