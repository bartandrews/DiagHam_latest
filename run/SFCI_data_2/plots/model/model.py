import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.patches import Circle
import matplotlib.patches as patches

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


if __name__ == "__main__":

    fig = plt.figure(figsize=(5, 5))
    ax = plt.gca()

    # grid
    for x in range(9):
        for y in range(9):
            if x == 4 and y == 4:
                clr = 'k'
            elif x == 4 or y == 4:
                idx = np.abs(4-x) if np.abs(x) != 4 else np.abs(4-y)
                clr = f'C{idx}'
            else:
                clr = 'gray'
            ax.plot(x-4, y-4, marker=".", markersize=10, c=clr)

    # circles
    radii = [1, np.sqrt(2), 2, np.sqrt(5), np.sqrt(8), 3, np.sqrt(10), np.sqrt(13), 4]
    for i in radii:
        if i in [1, 2, 3, 4]:
            clr = f'C{i}'
            lnst = '-'
        else:
            clr = 'gray'
            lnst = '--'
        ax.add_patch(Circle((0, 0), i, facecolor=None, edgecolor=clr, fill=None, zorder=-1, linewidth=0.5,
                            linestyle=lnst))

    # labels
    ax.text(0.2, 1.05, "1NN", c='C1')
    ax.text(0.2, 2.05, "3NN", c='C2')
    ax.text(0.2, 3.05, "6NN", c='C3')
    ax.text(0.2, 4.05, "9NN", c='C4')

    ax.text(0.05, 0.05, "$t_1=1$", c='C1')
    ax.text(1.4, -0.1, "$t_3$", c='C2')
    ax.text(2.4, -0.1, "$t_6$", c='C3')
    ax.text(3.4, -0.1, "$t_9$", c='C4')

    ax.text(-3.75, 3.35, "$\odot\;\mathbf{B}$", c='k', fontsize=14)

    # arrows
    style = "Simple, tail_width=0.5, head_width=4, head_length=4"
    kw = dict(arrowstyle=style)
    ax.add_patch(patches.FancyArrowPatch((0, 0), (1, 0), connectionstyle="arc3,rad=.25", zorder=0, color="C1", **kw))
    ax.add_patch(patches.FancyArrowPatch((0, 0), (2, 0), connectionstyle="arc3,rad=.25", zorder=0, color="C2", **kw))
    ax.add_patch(patches.FancyArrowPatch((0, 0), (3, 0), connectionstyle="arc3,rad=.25", zorder=0, color="C3", **kw))
    ax.add_patch(patches.FancyArrowPatch((0, 0), (4, 0), connectionstyle="arc3,rad=.25", zorder=0, color="C4", **kw))

    ax.add_patch(patches.FancyArrowPatch((0, 0), (0, 1), connectionstyle="arc3,rad=-0.25", zorder=0, color="C1", **kw))
    ax.add_patch(patches.FancyArrowPatch((0, 0), (0, 2), connectionstyle="arc3,rad=-0.25", zorder=0, color="C2", **kw))
    ax.add_patch(patches.FancyArrowPatch((0, 0), (0, 3), connectionstyle="arc3,rad=-0.25", zorder=0, color="C3", **kw))
    ax.add_patch(patches.FancyArrowPatch((0, 0), (0, 4), connectionstyle="arc3,rad=-0.25", zorder=0, color="C4", **kw))

    # axes
    ax.set_xlabel("$x/a$")
    ax.set_ylabel("$y/a$")

    plt.savefig(f"/home/bart/DiagHam_latest/run/SFCI_data_2/plots/model/model.png", bbox_inches='tight', dpi=300)
    plt.show()
