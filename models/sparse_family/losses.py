import torch
import torch.nn.functional as F

def family_loss(family_logits: torch.Tensor, y_family: torch.Tensor) -> torch.Tensor:
    """
    Cross-entropy loss for FAMILY prediction, ignoring samples with label -1.

    Args:
        family_logits: Tensor of shape (batch_size, num_family)
        y_family: Tensor of shape (batch_size,) with family indices or -1 for unknowns

    Returns:
        torch.Tensor: Scalar loss value
    """
    mask = (y_family != -1)

    if mask.any():
        return F.cross_entropy(family_logits[mask], y_family[mask])
    else:
        return torch.tensor(0.0, device=family_logits.device)
