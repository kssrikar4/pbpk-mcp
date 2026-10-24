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

def plot_multicurve_comparison(scenarios, output_path):
    if not isinstance(scenarios, (list, dict)):
        raise TypeError("scenarios must be a list or dict of simulation objects.")
    if not isinstance(output_path, str) or len(output_path.strip()) == 0:
        raise ValueError("Invalid output_path.")

    if isinstance(scenarios, dict):
        items = [{"label": k, "data": v} for k, v in scenarios.items()]
    else:
        items = scenarios

    palette = ["#0F172A", "#2563EB", "#D97706", "#DC2626", "#059669", "#7C3AED"]

    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=300)

    for idx, item in enumerate(items):
        lbl = item.get("label", "Scenario " + str(idx + 1))
        d = item.get("data", item)
        time_arr = np.array(d["time_h"], dtype=np.float64)
        c_arr = np.array(d["venous_plasma_mg_l"], dtype=np.float64)
        color = palette[idx % len(palette)]
        ax.plot(time_arr, c_arr, label=lbl, color=color, linewidth=1.75)

    ax.set_xlabel("Time (h)", fontsize=10, fontweight="medium", color="#111827")
    ax.set_ylabel("Plasma Concentration (mg/L)", fontsize=10, fontweight="medium", color="#111827")

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

def plot_auc_barchart(labels, aucs, output_path):
    if not isinstance(labels, list) or not isinstance(aucs, list):
        raise TypeError("labels and aucs must be lists.")
    if len(labels) != len(aucs) or len(labels) == 0:
        raise ValueError("labels and aucs must have matching non-zero lengths.")

    palette = ["#0F172A", "#2563EB", "#D97706", "#DC2626", "#059669"]

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=300)

    x_indices = np.arange(len(labels))
    bars = ax.bar(x_indices, aucs, color=[palette[i % len(palette)] for i in range(len(labels))], width=0.55, edgecolor="none")

    for bar, val in zip(bars, aucs):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, h + max(aucs) * 0.02, "{:.1f}".format(val), ha="center", va="bottom", fontsize=8.5, color="#111827", fontweight="medium")

    formatted_labels = [str(l).replace(" (", "\n(") for l in labels]
    ax.set_xticks(x_indices)
    ax.set_xticklabels(formatted_labels, fontsize=8.5, color="#111827", rotation=15, ha="right")
    ax.set_ylabel("AUC (mg·h/L)", fontsize=10, fontweight="medium", color="#111827")
    ax.set_ylim(0, max(aucs) * 1.15)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B7280")
    ax.spines["bottom"].set_color("#6B7280")
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.tick_params(axis="both", which="major", labelsize=9, colors="#374151", direction="out", length=3.5, width=0.8)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    return output_path

def plot_renal_function_vs_exposure(renal_metrics, exposures, output_path, xlabel="eGFR (mL/min)", ylabel="AUC (mg·h/L)"):
    if not isinstance(renal_metrics, list) or not isinstance(exposures, list):
        raise TypeError("renal_metrics and exposures must be lists.")
    if len(renal_metrics) != len(exposures) or len(renal_metrics) == 0:
        raise ValueError("renal_metrics and exposures must have matching non-zero lengths.")

    r_arr = np.array(renal_metrics, dtype=np.float64)
    e_arr = np.array(exposures, dtype=np.float64)

    fig, ax = plt.subplots(figsize=(6.5, 4.0), dpi=300)

    ax.scatter(r_arr, e_arr, color="#0F172A", s=40, zorder=3)

    sort_idx = np.argsort(r_arr)
    ax.plot(r_arr[sort_idx], e_arr[sort_idx], color="#2563EB", linewidth=1.75, linestyle="--", zorder=2)

    ax.set_xlabel(xlabel, fontsize=10, fontweight="medium", color="#111827")
    ax.set_ylabel(ylabel, fontsize=10, fontweight="medium", color="#111827")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B7280")
    ax.spines["bottom"].set_color("#6B7280")
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)

    ax.tick_params(axis="both", which="major", labelsize=9, colors="#374151", direction="out", length=3.5, width=0.8)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    return output_path
