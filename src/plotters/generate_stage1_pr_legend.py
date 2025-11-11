import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Reuse styling from the PR curve script to ensure consistency
from generate_stage1_pr_curve import (
    assembler_base_colors,
    model_shade_factors,
    BASELINE_SHADE_FACTOR,
    BASELINE_LINESTYLE,
    shade_hex_color,
    STAGE1_PR_PLOT_FOLDER,
)


def generate_legend_for_tools(tool_names: list[str], output_filename: str) -> str:
    """Create a standalone legend PDF for a given set of tools.

    tool_names should be a list like ["stringtie", "isoquant"] or ["stringtie", "scallop2"].
    Colors are taken from assembler_base_colors using:
      stringtie -> assembler1, isoquant/scallop2 -> assembler2
    """
    os.makedirs(STAGE1_PR_PLOT_FOLDER, exist_ok=True)
    out_path = os.path.join(STAGE1_PR_PLOT_FOLDER, output_filename)

    tool_to_assembler = {
        "stringtie": "assembler1",
        "isoquant": "assembler2",
        "scallop2": "assembler2",
    }
    pretty_name = {
        "stringtie": "StringTie",
        "isoquant": "IsoQuant",
        "scallop2": "Scallop2",
    }

    model_map = {
        "xgboost": "Telos-XGB",
        "randomforest": "Telos-RF"
    }
    # Collect handles per tool so we can place each tool on its own row
    tool_handles = {}
    model_order = ["xgboost", "randomforest"]

    for tool in tool_names:
        assembler_key = tool_to_assembler.get(tool, "assembler2")
        base_color = assembler_base_colors.get(assembler_key, "#333333")

        # Baseline entry
        baseline_color = shade_hex_color(base_color, BASELINE_SHADE_FACTOR)
        handles = [
            Line2D(
                [0], [0], color=baseline_color, linestyle=BASELINE_LINESTYLE, linewidth=1,
                label=f"{pretty_name.get(tool, tool.title())} Baseline",
            )
        ]

        # Model entries
        for model in model_order:
            shade = model_shade_factors.get(model, 1.0)
            model_color = shade_hex_color(base_color, shade)
            handles.append(
                Line2D(
                    [0], [0], color=model_color, linestyle="-", linewidth=2,
                    label=f"{pretty_name.get(tool, tool.title())} {model_map.get(model, model.title())}",
                )
            )

        tool_handles[tool] = handles

    # Create figure sized to the legend only (slightly taller to fit two rows)
    fig, ax = plt.subplots(figsize=(10, 1.6 if len(tool_names) <= 1 else 1.9))
    ax.axis("off")
    # Remove all outer margins
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Place each tool's legend in a separate row, no frame, compact spacing
    row_positions = [0.55, 0.45]
    legend_objs = []
    for idx, tool in enumerate(tool_names):
        handles = tool_handles.get(tool, [])
        if not handles:
            continue
        leg = ax.legend(
            handles=handles,
            loc="center",
            bbox_to_anchor=(0.5, row_positions[idx] if idx < len(row_positions) else 0.5),
            ncol=len(handles),
            frameon=False,
            fontsize=10,
            columnspacing=0.8,
            handlelength=2.0,
            handletextpad=0.4,
            borderaxespad=0.0,
        )
        ax.add_artist(leg)
        legend_objs.append(leg)

    # Render, compute tight bbox around legend objects only, and save with no padding
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    from matplotlib.transforms import Bbox
    if legend_objs:
        bboxes = [leg.get_window_extent(renderer) for leg in legend_objs]
        bbox = Bbox.union(bboxes)
        bbox_inches = bbox.transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(out_path, dpi=300, format="pdf", bbox_inches=bbox_inches, pad_inches=0)
    else:
        fig.savefig(out_path, dpi=300, format="pdf", bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    long_reads_path = generate_legend_for_tools(["stringtie", "isoquant"], "stage1_pr_legend_long_reads.pdf")
    print(f"Legend saved: {long_reads_path}")
    short_reads_path = generate_legend_for_tools(["stringtie", "scallop2"], "stage1_pr_legend_short_reads.pdf")
    print(f"Legend saved: {short_reads_path}")


