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
                 num_classes: int = 1000,
                 dropout: float = 0.3):
        super().__init__()

        self.fc1 = nn.Linear(input_size, hidden_dim, bias=False)
        self.dropout = nn.Dropout(p=dropout)
        self.output_layer = nn.Linear(hidden_dim, num_classes)

    def forward(self, x_sparse: torch.Tensor) -> torch.Tensor:
        """
        Forward pass with sparse input tensor.
        Args:
            x_sparse (torch.sparse.Tensor): Input matrix (batch_size x input_size)
        Returns:
            torch.Tensor: Genus logits (batch_size x num_classes)
        """
        x = torch.sparse.mm(x_sparse, self.fc1.weight.T)
        x = F.relu(x)
        x = self.dropout(x)
        return self.output_layer(x)
    
    

class SparseNN_clust(nn.Module):
    """
    Simple feedforward network for sparse VPF hit vectors.
    Predicts genus logits from a fixed-length sparse input vector.
    """
    def __init__(self,
                 input_size: int,
                 hidden_dim: int = 2024,
                 num_classes: int = 1000,
                 dropout: float = 0.3):
        super().__init__()
        self.fc1 = nn.Linear(input_size, hidden_dim, bias=False)
        self.dropout = nn.Dropout(p=dropout)
        self.output_layer = nn.Linear(hidden_dim, num_classes)

    def forward(self, x_sparse: torch.Tensor, return_embeddings: bool=False) -> torch.Tensor:
        """
        Forward pass with sparse input tensor.
        Args:
            x_sparse (torch.sparse.Tensor): Input matrix (batch_size x input_size)
        Returns:
            torch.Tensor: Genus logits (batch_size x num_classes)
        """
        x = torch.sparse.mm(x_sparse, self.fc1.weight.T)
        x = F.relu(x)
        x = self.dropout(x)
        if return_embeddings == True:
            return x
        return self.output_layer(x)
    
    

class ProtoNN(nn.Module):
    """
    
    """
    def __init__(self,
                 input_size: int,
                 hidden_dim: int = 2024,
                 num_classes: int = 1000,
                 dropout: float = 0.3,
                 normalize_embeddigns: bool = True):
        super().__init__()

        self.fc1 = nn.Linear(input_size, hidden_dim, bias=False)
        self.dropout = nn.Dropout(p=dropout)
        self.prototypes = nn.Parameter(torch.randn(num_classes, hidden_dim))
        self.normalize_embeddigns = normalize_embeddigns

        # Optional: initialize prototypes with Xavier normal for stability
        nn.init.xavier_normal_(self.prototypes)

    def forward(self, x_sparse: torch.Tensor,
                return_embedding: bool = False) -> torch.Tensor:
        
        x = torch.sparse.mm(x_sparse, self.fc1.weight.T)
        x = F.relu(x)
        x = self.dropout(x)

        if self.normalize_embeddigns:
            x = F.normalize(x, dim=1)
            prototypes = F.normalize(self.prototypes, dim=1)
        else:
            prototypes = self.prototypes

        if return_embedding:
            return x
        
        dists = torch.cdist(x, prototypes, p=2)**2
        logits = -dists
        
        return logits





