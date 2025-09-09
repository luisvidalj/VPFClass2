import torch
from torch.utils.data import DataLoader
from tqdm import tqdm
from .losses import genus_loss


def train_one_epoch(model: torch.nn.Module,
                    dataloader: DataLoader,
                    optimizer: torch.optim.Optimizer,
                    device: torch.device,
                    epoch: int = 0,
                    print_every: int = 10) -> float:
    """
    Trains the model for one epoch using the genus loss.

    Args:
        model: PyTorch model
        dataloader: DataLoader with collate_fn_genus
        optimizer: Optimizer instance
        device: CPU or CUDA
        epoch: Epoch index (for printing/logging)
        print_every: Print frequency (in batches)

    Returns:
        Average loss across the epoch
    """
    model.train()
    running_loss = 0.0
    correct, total = 0, 0

    loop = tqdm(enumerate(dataloader), total=len(dataloader), desc=f"Epoch {epoch}")

    for i, (x_batch, y_batch, _) in loop:
        x_batch = x_batch.to(device)
        y_batch = y_batch.to(device)

        optimizer.zero_grad()
        outputs = model(x_batch)
        loss = genus_loss(outputs, y_batch)

        loss.backward()
        optimizer.step()

        running_loss += loss.item()

        # Accuracy (solo donde y != -1)
        with torch.no_grad():
            mask = (y_batch != -1)
            preds = outputs[mask].argmax(dim=1)
            correct += (preds == y_batch[mask]).sum().item()
            total += mask.sum().item()

        if print_every and i % print_every == 0:
            acc = 100 * correct / total if total > 0 else 0
            loop.set_postfix(loss=loss.item(), acc=f"{acc:.2f}%")

    epoch_loss = running_loss / len(dataloader)
    print(f"Epoch {epoch} complete. Loss: {epoch_loss:.4f}. Accuracy: {100 * correct / total:.2f}%")

    return epoch_loss
