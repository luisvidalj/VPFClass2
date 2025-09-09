import torch

def evaluate_genus_accuracy(model, dataloader, device='cpu'):
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
    print(f"Test Accuracy: {accuracy:.4f}")
    return accuracy