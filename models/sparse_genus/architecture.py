import torch
import torch.nn as nn
import torch.nn.functional as F


class SparseNN(nn.Module):
    """
    Simple feedforward network for sparse VPF hit vectors.
    Predicts genus logits from a fixed-length sparse input vector.
    """
    def __init__(self,
                 input_size: int,
                 hidden_dim: int = 2024,
                 num_genus: int = 1000,
                 dropout: float = 0.3):
        super().__init__()

        self.fc1 = nn.Linear(input_size, hidden_dim, bias=False)
        self.dropout = nn.Dropout(p=dropout)
        self.output_layer = nn.Linear(hidden_dim, num_genus)

    def forward(self, x_sparse: torch.Tensor) -> torch.Tensor:
        """
        Forward pass with sparse input tensor.
        Args:
            x_sparse (torch.sparse.Tensor): Input matrix (batch_size x input_size)
        Returns:
            torch.Tensor: Genus logits (batch_size x num_genus)
        """
        x = torch.sparse.mm(x_sparse, self.fc1.weight.T)
        x = F.relu(x)
        x = self.dropout(x)
        return self.output_layer(x)
    



