import matplotlib.pyplot as plt


def robustness_heatmaps(square, smooth, path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)
    for ax, grid, title in zip(axes, (square, smooth), ("Square", "Smooth")):
        image = ax.imshow(grid, origin="lower", aspect="auto",
                          extent=[-10, 10, -10, 10], vmin=0.5, vmax=1)
        ax.set_title(title)
        ax.set_xlabel("Pulse duration error (%)")
    axes[0].set_ylabel(r"Detuning / $\Omega$ (%)")
    fig.colorbar(image, ax=axes, label="Noiseless average fidelity")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig
