
import torch
import torch.nn.functional as F
import pandas as pd

def collect_predictions_genus(
    model,
    dataloader,
    idx_to_genus,
    device='cuda',
    max_print=None,
    topk=3
):
    """
    Collects genus predictions and top-k confidence scores for later analysis.

    Args:
        model: Trained model
        dataloader: DataLoader with GenusDataset
        idx_to_genus: Mapping from index to genus name
        device: 'cpu' or 'cuda'
        max_print: Maximum number of samples to process (for display)
        topk: Number of top predictions to collect

    Returns:
        count_genus_none: How many samples had missing true genus
        count: Number of samples processed
        df_preds: DataFrame with columns ['Accession', 'real_genus', 'pred_genus', 'topk_genus', 'vpfclass_score']
    """
    model.eval()
    model.to(device)
    count = 0
    count_genus_none = 0
    records = []

    with torch.no_grad():
        for x_batch, y_batch, accessions in dataloader:
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            genus_logits = model(x_batch)
            genus_probs = F.softmax(genus_logits, dim=1)
            genus_preds = genus_logits.argmax(dim=1)

            for i in range(len(x_batch)):
                real_genus = idx_to_genus.get(y_batch[i].item(), 'None')
                pred_genus = idx_to_genus.get(genus_preds[i].item(), 'None')
                accession  = accessions[i]

                topk_gen_probs, topk_gen_idx = torch.topk(genus_probs[i], topk)
                topk_gen_info = [
                    f"{idx_to_genus.get(topk_gen_idx[j].item(), 'None')} ({topk_gen_probs[j].item():.2f})"
                    for j in range(topk)
                ]
                topk_gen_str = "; ".join(topk_gen_info)

                vpfclass_score = topk_gen_probs[0].item()  # Probabilidad de la predicción top-1

                if real_genus == 'None':
                    count_genus_none += 1

                print(f"Real: {real_genus:20s} | {accession}")
                print(f"Pred: {pred_genus:20s}")
                print(f"Top-{topk} Genus: {topk_gen_str}")
                print("-" * 60)

                records.append({
                    'Accession':      accession,
                    'real_genus':     real_genus,
                    'pred_genus':     pred_genus,
                    'topk_genus':     topk_gen_str,
                    'vpfclass_score': vpfclass_score
                })

                count += 1
                if max_print and count >= max_print:
                    return count_genus_none, count, pd.DataFrame(records)

    return count_genus_none, count, pd.DataFrame(records)




import torch
import torch.nn.functional as F
import pandas as pd
from sparse_genus.evaluate import compute_energy  # o del módulo donde la tengas


def collect_predictions_genus_e(
    model,
    dataloader,
    idx_to_genus,
    device='cuda',
    max_print=None,
    topk=3,
    apply_energy_filter=False,
    tau_energy=None,
    temperature=1.0
):
    """
    Collects genus predictions, top-k confidence scores, and optionally applies
    an energy-based filter to label uncertain samples as 'unknown_genus'.

    Args:
        model: Trained model
        dataloader: DataLoader with GenusDataset
        idx_to_genus: Mapping from index to genus name
        device: 'cpu' or 'cuda'
        max_print: Maximum number of samples to process (for display)
        topk: Number of top predictions to collect
        apply_energy_filter: If True, use energy score to assign 'unknown_genus'
        tau_energy: Energy threshold for unknown detection
        temperature: Temperature used in energy computation

    Returns:
        count_genus_none: Number of samples with missing true genus
        count: Total number of processed samples
        df_preds: DataFrame with columns:
            ['Accession', 'real_genus', 'pred_genus', 'topk_genus',
             'vpfclass_score', 'energy', 'is_unknown']
    """
    model.eval()
    model.to(device)
    count = 0
    count_genus_none = 0
    records = []

    with torch.no_grad():
        for x_batch, y_batch, accessions in dataloader:
            x_batch = x_batch.to(device)
            y_batch = y_batch.to(device)

            genus_logits = model(x_batch)
            genus_probs = F.softmax(genus_logits, dim=1)
            genus_preds = genus_logits.argmax(dim=1)

            # --- Energy computation ---
            if apply_energy_filter:
                energy_scores = compute_energy(genus_logits, temperature=temperature)
            else:
                energy_scores = torch.zeros(len(x_batch), device=device)

            for i in range(len(x_batch)):
                real_genus = idx_to_genus.get(y_batch[i].item(), 'None')
                pred_genus = idx_to_genus.get(genus_preds[i].item(), 'None')
                accession = accessions[i]

                # Energy & unknown flag
                energy = energy_scores[i].item()
                is_unknown = (
                    apply_energy_filter and (tau_energy is not None) and (energy > tau_energy)
                )

                if is_unknown:
                    pred_genus = 'unknown_genus'

                topk_gen_probs, topk_gen_idx = torch.topk(genus_probs[i], topk)
                topk_gen_info = [
                    f"{idx_to_genus.get(topk_gen_idx[j].item(), 'None')} ({topk_gen_probs[j].item():.2f})"
                    for j in range(topk)
                ]
                topk_gen_str = "; ".join(topk_gen_info)

                vpfclass_score = topk_gen_probs[0].item()  # Probabilidad top-1

                if real_genus == 'None':
                    count_genus_none += 1

                # --- Print (optional) ---
                print(f"Real: {real_genus:20s} | {accession}")
                if is_unknown:
                    print(f"Pred: {'unknown_genus':20s}  [Energy={energy:.2f} > {tau_energy:.2f}]")
                else:
                    print(f"Pred: {pred_genus:20s}  [Energy={energy:.2f}]")
                print(f"Top-{topk} Genus: {topk_gen_str}")
                print("-" * 60)

                records.append({
                    'Accession': accession,
                    'real_genus': real_genus,
                    'pred_genus': pred_genus,
                    'topk_genus': topk_gen_str,
                    'vpfclass_score': vpfclass_score,
                    'energy': energy,
                    'is_unknown': is_unknown
                })

                count += 1
                if max_print and count >= max_print:
                    return count_genus_none, count, pd.DataFrame(records)

    return count_genus_none, count, pd.DataFrame(records)


import numpy as np

@torch.no_grad()
def extract_embeddings(model, dataloader, device='cuda'):
    model.eval(); model.to(device)
    embs, y_true, accs, maxprob, preds = [], [], [], [], []

    for x_batch, y_batch, accessions in dataloader:
        x = x_batch.to(device)
        logits = model(x)
        p = F.softmax(logits, dim=1)
        mprob, pred = p.max(dim=1)
        h = model(x, return_embeddings= True)

        embs.append(h.cpu().numpy())
        y_true.append(y_batch.cpu().numpy())
        accs.extend(accessions)
        maxprob.append(mprob.cpu().numpy())
        preds.append(pred.cpu().numpy())

    E = np.vstack(embs)
    y = np.concatenate(y_true)
    A = np.array(accs, dtype=object)
    P = np.concatenate(maxprob)
    Yhat = np.concatenate(preds)
    return E, y, A, P, Yhat


def mask_confident(y, Yhat, P, tau_conf=0.8, require_correct=False):
    m = P >= tau_conf
    if require_correct:
        m = m & (y == Yhat) & (y != -1)
    else:
        m = m & (y != -1)
    return m


from sklearn.decomposition import PCA

def fit_pca(E, n_components=128, random_state=0):
    pca = PCA(n_components=n_components, random_state=random_state, svd_solver='auto')
    Z = pca.fit_transform(E)
    return pca, Z



from sklearn.mixture import GaussianMixture

def fit_gmm_per_class(Z, y, num_classes, 
                      kmax_by_size=lambda n: max(1, min(10, n // 30)),
                      cov_type='diag', random_state=0):
    """
    Z: embeddings (posiblemente PCA) de los ejemplos confiables
    y: labels correspondientes
    num_classes: total de clases
    kmax_by_size: función que define k_max en función de n_c
    """
    models = {}
    stats = {}  # para thresholds por componente
    for c in range(num_classes):
        idx = np.where(y == c)[0]
        if len(idx) == 0:
            continue
        Zc = Z[idx]
        n_c = len(Zc)

        if n_c == 1: # Que metemos como varianza (solo algo dif de 0)
            mean = Zc[0]
            cov = np.ones(Zc.shape[1]) * 1e-5
            gmm_stub = type("DummyGMM", (), {})()  # objeto vacío tipo GMM
            gmm_stub.n_components = 1
            gmm_stub.means_ = np.expand_dims(mean, axis=0)
            gmm_stub.covariances_ = np.expand_dims(cov, axis=0)
            models[c] = gmm_stub
            stats[c] = [{"q95": 0.0, "n": 1}]
            continue

        kmax = kmax_by_size(len(idx))
        best_bic, best_gmm = np.inf, None
        for k in range(1, kmax+1):
            gmm = GaussianMixture(n_components=k, covariance_type=cov_type, 
                                  reg_covar=1e-5, random_state=random_state)
            gmm.fit(Zc)
            bic = gmm.bic(Zc)
            if bic < best_bic:
                best_bic, best_gmm = bic, gmm
        models[c] = best_gmm

        # Distancias (Mahalanobis) intra-cluster para fijar umbrales por componente
        # Para diag-cov podemos calcular distancia rápida:
        means = best_gmm.means_
        covs  = best_gmm.covariances_  # (k, d) si diag
        # asignaciones
        resp = best_gmm.predict(Zc)
        comp_stats = []
        for j in range(best_gmm.n_components):
            Zcj = Zc[resp == j]
            if len(Zcj) == 0:
                comp_stats.append({'q95': np.inf, 'n': 0})
                continue
            diff = Zcj - means[j]
            inv_var = 1.0 / covs[j]
            # Mahalanobis^2 (diag): sum((diff^2) * inv_var, axis=1)
            d2 = np.sum((diff**2) * inv_var, axis=1)
            q95 = np.percentile(d2, 95)  # umbral dentro de comp
            comp_stats.append({'q95': float(q95), 'n': int(len(Zcj))})
        stats[c] = comp_stats
    return models, stats


import numpy as np

def mahalanobis2_diag(x, mean, diag_cov):
    inv_var = 1.0 / diag_cov
    diff = x - mean
    return float(np.sum((diff**2) * inv_var))

def route_with_clusters(model, dataloader, idx_to_genus,
                        pca, gmm_models, comp_stats,
                        tau_conf=0.8, device='cuda'):
    """
    Hybrid inference:
      - If max_softmax >= tau_conf -> use SparseNN prediction (route='sparse')
      - Otherwise:
          * compute min Mahalanobis² to any cluster (c,j)
          * if min_d2 <= q95_{c,j} -> assign that class (route='cluster')
          * else -> route='unknown'
    Also records the nearest cluster for every sample (cluster_pred, cluster_d2),
    even if the route is 'sparse'.
    """
    model.eval()
    model.to(device)
    records = []

    with torch.no_grad():
        for x_batch, y_batch, accessions in dataloader:
            # --- Forward pass ---
            x = x_batch.to(device)
            logits = model(x)
            probs = F.softmax(logits, dim=1)
            mprob, yhat = probs.max(dim=1)
            H = model(x, return_embeddings=True).cpu().numpy()  # embeddings
            Z = pca.transform(H) if pca is not None else H     # optional PCA

            # --- Loop over samples ---
            for i in range(len(x_batch)):
                acc = accessions[i]
                real = idx_to_genus.get(int(y_batch[i]), 'None')
                mp = float(mprob[i].cpu().item())
                pred = int(yhat[i].cpu().item())
                z = Z[i]

                # 🔹 Buscar SIEMPRE el cluster más cercano
                best_c, best_j, d2_min = None, None, np.inf
                for c, gmm in gmm_models.items():
                    if gmm is None:
                        continue
                    for j in range(gmm.n_components):
                        d2 = mahalanobis2_diag(z, gmm.means_[j], gmm.covariances_[j])
                        if d2 < d2_min:
                            d2_min, best_c, best_j = d2, c, j

                cluster_pred = idx_to_genus.get(best_c, 'None') if best_c is not None else None

                # 🔹 Decidir la ruta
                if mp >= tau_conf:
                    final_cls = pred
                    routed = 'sparse'
                else:
                    if best_c is not None:
                        q95 = comp_stats[best_c][best_j]['q95']
                        if d2_min <= q95:
                            final_cls = best_c
                            routed = 'cluster'
                        else:
                            final_cls = None
                            routed = 'unknown'
                    else:
                        final_cls = None
                        routed = 'unknown'

                final_genus = (
                    idx_to_genus.get(final_cls, 'unknown_genus')
                    if final_cls is not None else 'unknown_genus'
                )

                records.append({
                    'Accession': acc,
                    'real_genus': real,
                    'sparse_pred': idx_to_genus.get(pred, 'None'),
                    'sparse_conf': mp,
                    'route': routed,            # 'sparse' | 'cluster' | 'unknown'
                    'cluster_pred': cluster_pred,  # 👈 nuevo
                    'final_pred': final_genus,
                    'cluster_d2': d2_min           # 👈 nuevo (distancia mínima)
                })

    return pd.DataFrame(records)
