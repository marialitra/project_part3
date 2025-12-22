from libraries import torch, os, np, esm, SeqIO
import argparse

# -----------------------------
# 0. Configuration
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description='Generate protein embeddings using ESM2')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output .dat file')
    parser.add_argument('-model', '--model', type=str, default='esm2_t6_8M_UR50D', help='ESM2 model name')
    parser.add_argument('--batch-size', type=int, default=64, help='Batch size for processing')
    parser.add_argument('--max-len', type=int, default=1022, help='Maximum sequence length')
    
    args = parser.parse_args()
    
    FASTA = args.input
    OUTPUT_FILE = args.output
    BATCH_SIZE = args.batch_size
    MAX_LEN = args.max_len
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    MAX_ID_LEN = 50  # max characters per ID
    EMBED_DIM = 320  # esm2_t6_8M_UR50D

    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    # -----------------------------
    # 1. Load ESM model
    # -----------------------------
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model = model.to(DEVICE)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    # -----------------------------
    # 2. Load sequences
    # -----------------------------
    data = []
    for record in SeqIO.parse(FASTA, "fasta"):
        # For protein's id we want from >sp|A7IC07|RBFA_XANP2
        # to take the middle part (A7IC07)
        header_parts = record.id.split('|')
        if len(header_parts) == 3:
            protein_id = header_parts[1]
        else:
            protein_id = record.id  # fallback if header is not standard

        seq = str(record.seq)[:MAX_LEN]
        data.append((protein_id, seq))

    data.sort(key=lambda x: len(x[1]), reverse=True)
    N = len(data)

    # -----------------------------
    # 3. Prepare memmap for embeddings + IDs
    # -----------------------------
    # Total columns: EMBED_DIM + MAX_ID_LEN (IDs stored as float32 bytes)
    total_cols = EMBED_DIM + MAX_ID_LEN
    embedding_map = np.memmap(
        OUTPUT_FILE,
        dtype='float32',
        mode='w+',
        shape=(N, total_cols)
    )

    # -----------------------------
    # 4. Inference and write to memmap
    # -----------------------------
    def id_to_float_array(protein_id, max_len=MAX_ID_LEN):
        """Convert string ID to float32 array of fixed length."""
        b = protein_id.encode('ascii', errors='replace')[:max_len]  # bytes
        arr = np.zeros(max_len, dtype='float32')
        arr[:len(b)] = np.frombuffer(b, dtype=np.uint8)
        return arr

    with torch.no_grad():
        idx = 0
        for i in range(0, N, BATCH_SIZE):
            batch = data[i:i+BATCH_SIZE]
            labels, strs, tokens = batch_converter(batch)
            tokens = tokens.to(DEVICE)

            results = model(tokens, repr_layers=[6])
            reps = results["representations"][6]

            for j, tok in enumerate(tokens):
                valid = (tok != alphabet.padding_idx) & \
                        (tok != alphabet.cls_idx) & \
                        (tok != alphabet.eos_idx)

                emb = reps[j][valid].mean(dim=0).cpu().numpy()
                id_arr = id_to_float_array(batch[j][0])

                # Store: [embedding | id_array]
                embedding_map[idx, :EMBED_DIM] = emb
                embedding_map[idx, EMBED_DIM:] = id_arr
                idx += 1

    embedding_map.flush()
    print(f"Saved embeddings+IDs to {OUTPUT_FILE}")
    print(f"Shape: {embedding_map.shape} (N x {total_cols})")

if __name__ == "__main__":
    main()
