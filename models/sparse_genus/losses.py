import torch
import torch.nn.functional as F

def genus_loss(genus_logits: torch.Tensor, y_genus: torch.Tensor) -> torch.Tensor:
    """
    Cross-entropy loss for genus prediction, ignoring samples with label -1.

    Args:
        genus_logits: Tensor of shape (batch_size, num_genus)
        y_genus: Tensor of shape (batch_size,) with genus indices or -1 for unknowns

    Returns:
        torch.Tensor: Scalar loss value
    """
    mask = (y_genus != -1)

    if mask.any():
        return F.cross_entropy(genus_logits[mask], y_genus[mask])
    else:
        return torch.tensor(0.0, device=genus_logits.device)
