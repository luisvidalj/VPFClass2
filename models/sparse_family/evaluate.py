import torch

def evaluate_family_accuracy(model, dataloader, device: str = "cpu") -> float:
    """
    Evaluate accuracy of a family classification model.

    Args:
        model: Trained PyTorch model
        dataloader: DataLoader yielding (x_batch, y_family, accession)
        device: 'cpu' or 'cuda'

    Returns:
        Accuracy (float, 0–1)
    """
    model.eval()
    model.to(device)

    correct = 0
    total = 0

    with torch.no_grad():
        for x_batch, y_batch, _ in dataloader:
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            logits = model(x_batch)
            preds = logits.argmax(dim=1)

            mask = (y_batch != -1)
            correct += (preds[mask] == y_batch[mask]).sum().item()
            total += mask.sum().item()

    accuracy = correct / total if total > 0 else 0.0
    print(f"[EVAL][FAMILY] Test Accuracy: {accuracy:.4f}")
    return accuracy
