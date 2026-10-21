import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def plot_pbpk_curves(sim_result, output_path, target_organ="Liver"):
    if not isinstance(sim_result, dict):
        raise TypeError("sim_result must be a dictionary.")
    if "time_h" not in sim_result or "venous_plasma_mg_l" not in sim_result:
        raise ValueError("sim_result missing required keys: 'time_h' and 'venous_plasma_mg_l'.")
    if not isinstance(output_path, str) or len(output_path.strip()) == 0:
        raise ValueError("Invalid output_path.")

    time_arr = np.array(sim_result["time_h"], dtype=np.float64)
    plasma_c = np.array(sim_result["venous_plasma_mg_l"], dtype=np.float64)

    organs_dict = sim_result.get("organ_concentrations_mg_l", {})
    target_c = None
    if target_organ in organs_dict:
        target_c = np.array(organs_dict[target_organ], dtype=np.float64)

    fig, ax = plt.subplots(figsize=(6.5, 4.0), dpi=300)

    color_plasma = "#111827"
    color_target = "#0D9488"

    ax.plot(time_arr, plasma_c, color=color_plasma, linewidth=1.75, label="Plasma")

    if target_c is not None:
        ax.plot(time_arr, target_c, color=color_target, linewidth=1.5, linestyle="--", label=target_organ)

    ax.set_xlabel("Time (h)", fontsize=10, fontweight="medium", color="#111827")
    ax.set_ylabel("Concentration (mg/L)", fontsize=10, fontweight="medium", color="#111827")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B7280")
    ax.spines["bottom"].set_color("#6B7280")
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.tick_params(axis="both", which="major", labelsize=9, colors="#374151", direction="out", length=3.5, width=0.8)

    ax.legend(frameon=False, fontsize=9, loc="upper right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    return output_path
