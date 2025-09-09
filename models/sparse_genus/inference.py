
import torch
import torch.nn.functional as F
import pandas as pd

def collect_predictions_genus(
    model,
    dataloader,
    idx_to_genus,
    device='cuda',
    max_print=None,
    topk=3
):
    """
    Collects genus predictions and top-k confidence scores for later analysis.

    Args:
        model: Trained model
        dataloader: DataLoader with GenusDataset
        idx_to_genus: Mapping from index to genus name
        device: 'cpu' or 'cuda'
        max_print: Maximum number of samples to process (for display)
        topk: Number of top predictions to collect

    Returns:
        count_genus_none: How many samples had missing true genus
        count: Number of samples processed
        df_preds: DataFrame with columns ['Accession', 'real_genus', 'pred_genus', 'topk_genus', 'vpfclass_score']
    """
    model.eval()
    model.to(device)
    count = 0
    count_genus_none = 0
    records = []

    with torch.no_grad():
        for x_batch, y_batch, accessions in dataloader:
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            genus_logits = model(x_batch)
            genus_probs = F.softmax(genus_logits, dim=1)
            genus_preds = genus_logits.argmax(dim=1)

            for i in range(len(x_batch)):
                real_genus = idx_to_genus.get(y_batch[i].item(), 'None')
                pred_genus = idx_to_genus.get(genus_preds[i].item(), 'None')
                accession  = accessions[i]

                topk_gen_probs, topk_gen_idx = torch.topk(genus_probs[i], topk)
                topk_gen_info = [
                    f"{idx_to_genus.get(topk_gen_idx[j].item(), 'None')} ({topk_gen_probs[j].item():.2f})"
                    for j in range(topk)
                ]
                topk_gen_str = "; ".join(topk_gen_info)

                vpfclass_score = topk_gen_probs[0].item()  # Probabilidad de la predicción top-1

                if real_genus == 'None':
                    count_genus_none += 1

                print(f"Real: {real_genus:20s} | {accession}")
                print(f"Pred: {pred_genus:20s}")
                print(f"Top-{topk} Genus: {topk_gen_str}")
                print("-" * 60)

                records.append({
                    'Accession':      accession,
                    'real_genus':     real_genus,
                    'pred_genus':     pred_genus,
                    'topk_genus':     topk_gen_str,
                    'vpfclass_score': vpfclass_score
                })

                count += 1
                if max_print and count >= max_print:
                    return count_genus_none, count, pd.DataFrame(records)

    return count_genus_none, count, pd.DataFrame(records)

