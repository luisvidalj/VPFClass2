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




def compute_energy(logits: torch.Tensor, temperature: float = 1.0) -> torch.Tensor:
    """
    Compute the energy score for each sample.
    Lower energy → more confident (in-distribution),
    Higher energy → likely OOD (unknown class).

    Args:
        logits: torch.Tensor of shape (batch_size, num_classes)
        temperature: scaling factor for logits (default=1.0)
    Returns:
        torch.Tensor of shape (batch_size,) with energy scores.
    """
    return -temperature * torch.logsumexp(logits / temperature, dim=1)


import torch
import torch.nn.functional as F

def evaluate_genus_accuracy_e(model,
                            dataloader,
                            device='cpu',
                            compute_energy_scores: bool = False,
                            temperature: float = 1.0,
                            tau_energy: float | None = None):
    """
    Evaluate model accuracy, optionally computing energy-based OOD detection.

    Args:
        model: trained model (SparseNN)
        dataloader: DataLoader (validation or test)
        device: 'cpu' or 'cuda'
        compute_energy_scores: if True, also compute energy and unknown flags
        temperature: scaling factor for energy score
        tau_energy: optional threshold to flag unknown samples
    Returns:
        accuracy (float)
        optionally also (energies, is_unknown_flags)
    """
    model.eval()
    model.to(device)

    correct, total = 0, 0
    all_energies, all_unknowns = [], []

    with torch.no_grad():
        for x_batch, y_batch, _ in dataloader:
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            logits = model(x_batch)
            preds = logits.argmax(dim=1)

            mask = (y_batch != -1)
            correct += (preds[mask] == y_batch[mask]).sum().item()
            total += mask.sum().item()

            # ---------- Energy-based OOD ----------
            if compute_energy_scores:
                energy = compute_energy(logits, temperature=temperature)
                all_energies.extend(energy.cpu().tolist())
                if tau_energy is not None:
                    all_unknowns.extend((energy > tau_energy).cpu().tolist())

    accuracy = correct / total if total > 0 else 0.0
    print(f"Test Accuracy: {accuracy:.4f}")

    if compute_energy_scores:
        if tau_energy is not None:
            n_unknown = sum(all_unknowns)
            print(f"Flagged as unknown: {n_unknown}/{len(all_unknowns)} samples "
                  f"({100*n_unknown/len(all_unknowns):.2f}%)")
        return accuracy, all_energies, all_unknowns if tau_energy is not None else all_energies

    return accuracy
