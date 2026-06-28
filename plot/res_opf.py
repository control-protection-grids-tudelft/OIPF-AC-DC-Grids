import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

# ============================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False
})

# ============================================================
res_branch_og = np.array([
    [1,2,1,1,0.99],
    [1,4,1,1,71.52],
    [1,5,0,0,0],
    [2,3,1,1,21.43],
    [2,4,0,0,0],
    [2,5,0,0,0],
    [3,4,0,0,0],
    [3,5,1,1,75.38],
    [4,5,1,1,6.71]
], dtype=float)

res_bus_og = np.array([
    [1,71.65],[2,20.37],[3,53.28],[4,79.10],[5,69.28]
], dtype=float)

res_branch_opf = np.array([
    [1,2,1,1,18.16],
    [1,4,0,1,0],
    [1,5,1,1,59.38],
    [2,3,1,1,18.45],
    [2,4,0,0,0],
    [2,5,1,1,85.14],
    [3,4,0,0,0],
    [3,5,0,1,0],
    [4,5,1,1,78.32]
], dtype=float)

res_bus_opf = np.array([
    [1,76.48],[2,47.15],[3,18.40],[4,79.10],[5,69.28]
], dtype=float)

res_branch_iopf = np.array([
    [1,2,1,1,44.47],
    [1,4,0,1,0],
    [1,5,0,0,0],
    [2,3,1,1,13.14],
    [2,4,1,1,96.87],
    [2,5,0,0,0],
    [3,4,0,0,0],
    [3,5,1,0,0],
    [4,5,1,1,20.96]
], dtype=float)

res_bus_iopf = np.array([
    [1,44.23],[2,35.87],[3,13.11],[4,79.10],[5,21.01]
], dtype=float)

# ============================================================
res_branch_comm = np.array([
    [1,2,1,1,65],
    [1,4,1,1,15],
    [2,3,1,1,0],
    [2,5,1,1,45],
    [3,5,1,1,0],
    [4,5,1,1,0]
], dtype=float)

res_bus_comm = np.array([
    [1,80.0],[2,20.0],[3,0.0],[4,15.0],[5,45.0]
], dtype=float)

# ============================================================
pos = {
    1:(0,2), 2:(0,1), 3:(0,0),
    4:(2.5,1.5), 5:(2.5,0.5)
}

# ============================================================
cmap = LinearSegmentedColormap.from_list(
    "soft_orange", [(0,"#FFF4EE"),(0.5,"#FD742D"),(1,"#C94A0F")]
)

cmap_comm = LinearSegmentedColormap.from_list(
    "soft_green", [(0, "#EEFDF5"), (0.5, "#185A56"), (1, "#0B3D3A")]
)

# ============================================================
def build_graph(branch, bus):
    G = nx.Graph()
    for r in bus:
        G.add_node(int(r[0]), pinj=r[1])
    for r in branch:
        G.add_edge(int(r[0]), int(r[1]),
                   sw_from=int(r[2]),
                   sw_to=int(r[3]),
                   power=r[4])
    return G


def draw_graph(ax, G, norm, pmin, prange, cmap_used,
               node_facecolor="#fdfaf6",
               node_edgecolor="#2b2b2b",
               switch_color="#6C47E5",
               open_edge_color="#9A9A9A",
               show_switch=True,
               show_flow_arrows=False,
               arrow_color="#C2A060",
               node_size_min=520,
               node_size_max=1500):
    # edges
    for u, v, d in G.edges(data=True):
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        closed = (d['sw_from'] == 1 and d['sw_to'] == 1)

        ax.plot([x1, x2], [y1, y2],
                linestyle='-' if closed else (0, (4, 3)),
                color=cmap_used(norm(d['power'])) if closed else open_edge_color,
                linewidth=1.5,
                zorder=1)

        # flow direction arrows
        if show_flow_arrows and abs(d['power']) > 1e-12:
            if d['power'] > 0:
                xs, ys = x1, y1
                xe, ye = x2, y2
            else:
                xs, ys = x2, y2
                xe, ye = x1, y1

            a1, a2 = 0.42, 0.58
            xa = xs + a1 * (xe - xs)
            ya = ys + a1 * (ye - ys)
            xb = xs + a2 * (xe - xs)
            yb = ys + a2 * (ye - ys)

            ax.annotate(
                "",
                xy=(xb, yb),
                xytext=(xa, ya),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=arrow_color,
                    lw=1.2,
                    mutation_scale=12,
                    shrinkA=0,
                    shrinkB=0
                ),
                zorder=6
            )

    # nodes
    nodes = list(G.nodes())
    pinj = np.array([G.nodes[n]['pinj'] for n in nodes], dtype=float)

    if prange < 1e-12:
        sizes = np.full_like(pinj, (node_size_min + node_size_max) / 2.0)
    else:
        sizes = node_size_min + (node_size_max - node_size_min) * (pinj - pmin) / prange

    nx.draw_networkx_nodes(
        G, pos, nodelist=nodes,
        node_size=sizes,
        node_color=node_facecolor,
        edgecolors=node_edgecolor,
        linewidths=2,
        ax=ax
    )

    nx.draw_networkx_labels(
        G, pos,
        labels={n: rf"$\#{n}$" for n in nodes},
        font_size=13,
        ax=ax
    )

    # switches
    if show_switch:
        for u, v, d in G.edges(data=True):
            x1, y1 = pos[u]
            x2, y2 = pos[v]

            for a, sw in [(0.22, d['sw_from']), (0.78, d['sw_to'])]:
                ax.scatter(
                    x1 + a * (x2 - x1),
                    y1 + a * (y2 - y1),
                    s=120, marker='*',
                    facecolors=switch_color if sw == 1 else 'white',
                    edgecolors=switch_color,
                    zorder=5
                )

    ax.set_aspect('auto')
    ax.axis('off')

    margin = 0.15
    ax.set_xlim(-0.45 - margin, 3.15 + margin)
    ax.set_ylim(-0.05 - margin, 2.1 + margin)

# ============================================================
Gs_power = [
    build_graph(res_branch_og, res_bus_og),
    build_graph(res_branch_opf, res_bus_opf),
    build_graph(res_branch_iopf, res_bus_iopf)
]

G_comm = build_graph(res_branch_comm, res_bus_comm)

# ============================================================
all_p = [d['power'] for G in Gs_power for _, _, d in G.edges(data=True)]
norm = Normalize(min(all_p), max(all_p), clip=True)

all_pinj = []
for G in Gs_power:
    all_pinj.extend([G.nodes[n]['pinj'] for n in G.nodes()])

global_pmin = float(np.min(all_pinj))
global_pmax = float(np.max(all_pinj))
global_prange = global_pmax - global_pmin

all_comm_p = [d['power'] for _, _, d in G_comm.edges(data=True)]
norm_comm = Normalize(min(all_comm_p), max(all_comm_p), clip=True)

all_comm_node = [G_comm.nodes[n]['pinj'] for n in G_comm.nodes()]
comm_pmin = float(np.min(all_comm_node))
comm_pmax = float(np.max(all_comm_node))
comm_prange = comm_pmax - comm_pmin

# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(25, 5.2))
plt.subplots_adjust(top=0.86, bottom=0.13, right=0.75, wspace=0.25)

titles = [
    r"$\mathrm{Original \ Elec.\ Net.}$",
    r"$\mathrm{Optimized\ Elec.\ Net.\ (OPF)}$",
    r"$\mathrm{Optimized\ Elec.\ Net.\ (OIPF)}$",
    r"$\mathrm{Comm.\ Net.\ (OIPF)}$"
]

shift_left = 0.03
for i in range(3):
    pos_i = axes[i].get_position()
    axes[i].set_position([
        pos_i.x0 - shift_left,
        pos_i.y0,
        pos_i.width,
        pos_i.height
    ])

# electric networks
for ax, G, t in zip(axes[:3], Gs_power, titles[:3]):
    draw_graph(
        ax, G, norm, global_pmin, global_prange, cmap_used=cmap,
        node_size_min=520, node_size_max=1500
    )
    ax.text(0.5, -0.015, t, transform=ax.transAxes,
            ha='center', va='bottom', fontsize=13)

# coommunication network
pos3 = axes[3].get_position()
axes[3].set_position([
    pos3.x0 + 0.08,
    pos3.y0,
    pos3.width,
    pos3.height
])

draw_graph(
    axes[3], G_comm, norm_comm, comm_pmin, comm_prange, cmap_used=cmap_comm,
    node_facecolor="#F2FBF7",
    node_edgecolor="#185A56",
    switch_color="#185A56",
    open_edge_color="#A8C8B8",
    show_switch=False,
    show_flow_arrows=True,
    arrow_color="#C2A060",
    node_size_min=8,
    node_size_max=1800
)

axes[3].text(0.5, -0.015, titles[3], transform=axes[3].transAxes,
             ha='center', va='bottom', fontsize=13)


cax = fig.add_axes([0.12, 0.87, 0.39, 0.028])
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
cbar.set_label(r"$\mathrm{Power\ Flow\ [MW]}$", fontsize=13)
cbar.ax.xaxis.set_label_position('top')
cbar.ax.xaxis.tick_top()
cbar.ax.tick_params(labelsize=12)

ticks = np.linspace(min(all_p), np.floor(max(all_p)), 5)
cbar.set_ticks(ticks)
cbar.set_ticklabels([rf"${x:g}$" for x in ticks])

cax_comm = fig.add_axes([0.88, 0.2, 0.01, 0.5])
sm_comm = plt.cm.ScalarMappable(cmap=cmap_comm, norm=norm_comm)
cbar_comm = fig.colorbar(sm_comm, cax=cax_comm, orientation='vertical')
cbar_comm.set_label(r"$\mathrm{Information\ Flow\ [MB/s]}$", fontsize=13)
cbar_comm.ax.tick_params(labelsize=12)

ticks_comm = np.linspace(min(all_comm_p), np.floor(max(all_comm_p)), 6)
cbar_comm.set_ticks(ticks_comm)
cbar_comm.set_ticklabels([rf"${x:g}$" for x in ticks_comm])

legend_power = [
    Line2D([0], [0], color="#C94A0F", lw=1.8,
           label=r"$\mathrm{Closed\ Branch}$"),
    Line2D([0], [0], color="#9A9A9A", lw=1.8, linestyle=(0, (4, 3)),
           label=r"$\mathrm{Open\ Branch}$"),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='#6C47E5',
           markeredgecolor='#6C47E5', markersize=10,
           label=r"$\mathrm{Closed\ Breaker}$"),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='white',
           markeredgecolor='#6C47E5', markersize=10,
           label=r"$\mathrm{Open\ Breaker}$"),
    Line2D([0], [0], marker='o', color='w',
           markerfacecolor='#fdfaf6',
           markeredgecolor='#2b2b2b',
           markersize=13,
           label=r"$\mathrm{DC\ Bus}$")
]

fig.legend(handles=legend_power, loc='lower left',
           ncol=5, frameon=False, fontsize=11.5,
           bbox_to_anchor=(0.08, 0.03))

legend_comm = [
    Line2D([0], [0],
           color="#185A56", lw=1.8,
           label=r"$\mathrm{Comm.\ Link}$"),

    Line2D([0], [0],
           marker='o', color='w',
           markerfacecolor='#F2FBF7',
           markeredgecolor='#185A56',
           markersize=13,
           label=r"$\mathrm{Comm.\ Node}$"),

    Line2D([0], [0],
           color="#C2A060", lw=1.2,
           marker='>',
           markerfacecolor="#C2A060",
           markeredgecolor="#C2A060",
           markersize=8,
           label=r"$\mathrm{Pos.\ Direction}$")
]

fig.legend(handles=legend_comm, loc='lower right',
           ncol=3, frameon=False, fontsize=11.5,
           bbox_to_anchor=(0.98, 0.03))

# ============================================================
plt.show()