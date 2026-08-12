import os, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import nnls
os.makedirs("figures", exist_ok=True); os.makedirs("results", exist_ok=True)
rng = np.random.default_rng(0)
cells = ["T CD8", "T CD4", "B cell", "NK", "Macrophage", "Neutrophil"]
G = 200                                          # marker genes
S = np.abs(rng.normal(1, 1, (G, len(cells))))    # signature matrix (genes x cell types)
for j in range(len(cells)):                      # each cell type over-expresses a marker block
    S[j*30:(j+1)*30, j] += 6
f_true = rng.dirichlet(np.ones(len(cells)))
bulk = S @ f_true + rng.normal(0, 0.3, G)        # simulated bulk expression
f_est, _ = nnls(S, bulk); f_est = f_est / f_est.sum()
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
x = np.arange(len(cells)); w = 0.38
ax[0].bar(x-w/2, f_true, w, label="true", color="#4C72B0")
ax[0].bar(x+w/2, f_est, w, label="estimated", color="#C44E52")
ax[0].set_xticks(x); ax[0].set_xticklabels(cells, rotation=25, fontsize=8); ax[0].set_ylabel("fraction"); ax[0].set_title("Deconvolved cell-type fractions"); ax[0].legend(fontsize=8)
ax[1].scatter(f_true, f_est, s=60, color="#8172B3"); lim=max(f_true.max(),f_est.max())*1.1
ax[1].plot([0,lim],[0,lim],"k--"); ax[1].set_xlabel("true fraction"); ax[1].set_ylabel("estimated fraction")
r = np.corrcoef(f_true, f_est)[0,1]; ax[1].set_title(f"Estimated vs true (r={r:.3f})")
fig.suptitle("Immune-cell deconvolution of bulk expression (demo data)"); fig.tight_layout(); fig.savefig("figures/demo.png", dpi=140)
open("results/summary.csv","w").write("cell_type,true,estimated\n"+"\n".join(f"{c},{t:.3f},{e:.3f}" for c,t,e in zip(cells,f_true,f_est))+"\n")
print(f"corr(true,est)={r:.3f}"); print("ok")
