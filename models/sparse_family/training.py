# sparse_family/training.py
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
from typing import Optional, Callable


def train_one_epoch(model: torch.nn.Module,
                    dataloader: DataLoader,
                    optimizer: torch.optim.Optimizer,
                    device: torch.device,
                    epoch: int = 0,
                    print_every: int = 10,
                    loss_fn: Optional[Callable[[torch.Tensor, torch.Tensor], torch.Tensor]] = None) -> float:
    """
    Train for one epoch (FAMILY task by defecto, pero agnóstico si pasas otro loss_fn).

    Args:
        model: PyTorch model
        dataloader: DataLoader con collate_fn_family (o equivalente)
        optimizer: Optimizer instance
        device: CPU o CUDA
        epoch: Índice de época (para logging)
        print_every: Frecuencia de impresión (en batches)
        loss_fn: Si None, usa CrossEntropyLoss(ignore_index=-1)

    Returns:
        Pérdida media de la época.
    """
    model.train()
    running_loss = 0.0
    correct, total = 0, 0

    if loss_fn is None:
        # Ignora ejemplos con etiqueta -1 (no asignables)
        loss_fn = torch.nn.CrossEntropyLoss(ignore_index=-1)

    loop = tqdm(enumerate(dataloader), total=len(dataloader), desc=f"Epoch {epoch}")

    for i, (x_batch, y_batch, _) in loop:
        x_batch = x_batch.to(device)
        y_batch = y_batch.to(device)

        optimizer.zero_grad()
        outputs = model(x_batch)  # logits (B, num_classes)
        loss = loss_fn(outputs, y_batch)

        loss.backward()
        optimizer.step()

        running_loss += loss.item()

        # Accuracy (solo donde y != -1)
        with torch.no_grad():
            mask = (y_batch != -1)
            if mask.any():
                preds = outputs[mask].argmax(dim=1)
                correct += (preds == y_batch[mask]).sum().item()
                total += mask.sum().item()

        if print_every and (i % print_every == 0):
            acc = 100.0 * correct / total if total > 0 else 0.0
            loop.set_postfix(loss=f"{loss.item():.4f}", acc=f"{acc:.2f}%")

    epoch_loss = running_loss / max(1, len(dataloader))
    acc_epoch = 100.0 * correct / total if total > 0 else 0.0
    print(f"Epoch {epoch} complete. Loss: {epoch_loss:.4f}. Accuracy: {acc_epoch:.2f}%")

    return epoch_loss
