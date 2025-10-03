
import torch
import torch.nn.functional as F
import pandas as pd

def collect_predictions_family(
    model,
    dataloader,
    idx_to_family,
    device: str = 'cuda',
    max_print: int | None = None,
    topk: int = 3
):
    """
    Collects FAMILY predictions and top-k confidence scores for later analysis.

    Args:
        model: Trained model
        dataloader: DataLoader with FamilyDataset
        idx_to_family: Mapping from index to family name
        device: 'cpu' or 'cuda'
        max_print: Maximum number of samples to print/process (for display)
        topk: Number of top predictions to collect

    Returns:
        count_family_none: How many samples had missing true family
        count: Number of samples processed
        df_preds: DataFrame with columns ['Accession', 'real_family', 'pred_family', 'topk_family', 'vpfclass_score']
    """
    model.eval()
    model.to(device)
    count = 0
    count_family_none = 0
    records = []

    with torch.no_grad():
        for x_batch, y_batch, accessions in dataloader:
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            family_logits = model(x_batch)
            family_probs = F.softmax(family_logits, dim=1)
            family_preds = family_logits.argmax(dim=1)

            for i in range(len(x_batch)):
                real_family = idx_to_family.get(y_batch[i].item(), 'None')
                pred_family = idx_to_family.get(family_preds[i].item(), 'None')
                accession   = accessions[i]

                topk_fam_probs, topk_fam_idx = torch.topk(family_probs[i], topk)
                topk_fam_info = [
                    f"{idx_to_family.get(topk_fam_idx[j].item(), 'None')} ({topk_fam_probs[j].item():.2f})"
                    for j in range(topk)
                ]
                topk_fam_str = "; ".join(topk_fam_info)

                vpfclass_score = topk_fam_probs[0].item()  # Probabilidad de la predicción top-1

                if real_family == 'None':
                    count_family_none += 1

                print(f"Real: {real_family:20s} | {accession}")
                print(f"Pred: {pred_family:20s}")
                print(f"Top-{topk} Family: {topk_fam_str}")
                print("-" * 60)

                records.append({
                    'Accession':      accession,
                    'real_family':    real_family,
                    'pred_family':    pred_family,
                    'topk_family':    topk_fam_str,
                    'vpfclass_score': vpfclass_score
                })

                count += 1
                if max_print and count >= max_print:
                    return count_family_none, count, pd.DataFrame(records)

    return count_family_none, count, pd.DataFrame(records)
