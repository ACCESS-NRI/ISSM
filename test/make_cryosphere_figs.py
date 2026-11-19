"""Figures for gpu_cryosphere_workshop.tex, written to figs_cryosphere/ as PDF.

  python make_cryosphere_figs.py [charts|map|all|breakdown|transient|mps|iterations]     (default all)

'charts' uses numbers from the AIS3 solver tests (ACCESS-AIS3/config/solver_*_test.py, outlogs in
ACCESS-AIS3/execution_gpu_test/, matrix lab logs in matrix_lab/). 'map' reads the AIS3 SSA model and
the dumped Picard-1 matrix, so it needs the pyissm environment (~3 min, ~6 GB).
Palette: dataviz reference categorical slots 1-2 (validated on white) and its blue sequential ramp.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figs_cryosphere')
BLUE, ORANGE = '#2a78d6', '#eb6834'            # categorical slots 1, 2
INK, INK2, MUTED = '#0b0b0b', '#52514e', '#898781'
GRID, AXIS = '#e1e0d9', '#c3c2b7'
SEQ = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5',
       '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']

plt.rcParams.update({
    # Drawn at the width they occupy on the slide (16:9 Beamer, ~5.5 in text width), so the
    # 9 pt text is not scaled down below the slide's own body text.
    'font.family': 'DejaVu Sans', 'font.size': 9, 'text.color': INK,
    'axes.edgecolor': AXIS, 'axes.linewidth': 0.8, 'axes.labelcolor': INK2,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
    'axes.spines.top': False, 'axes.spines.right': False,
    'legend.frameon': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03,
})


def grid(ax, axis='y'):
    ax.grid(axis=axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print('wrote', os.path.join(OUT, name))


def _breakdown(rows, xmax, step, name, figsize=(4.1, 1.5)):
    """Stacked bars: linear solves vs everything else, MUMPS (48 cores) against the GPU."""
    fig, ax = plt.subplots(figsize=figsize)
    for i, (lab, solve, rest, note) in enumerate(rows[::-1]):
        ax.barh(i, solve, height=0.42, color=BLUE, edgecolor='white', linewidth=1.5)
        ax.barh(i, rest, left=solve, height=0.42, color=ORANGE, edgecolor='white', linewidth=1.5)
        ax.text(solve + rest + xmax * 0.016, i, note, va='center', color=INK)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows[::-1]])
    ax.tick_params(axis='y', length=0)
    ax.set_xlim(0, xmax)
    ax.set_xticks(np.arange(0, xmax * 0.85, step))
    ax.set_xlabel('core solution time (s)')
    grid(ax, 'x')
    ax.spines['left'].set_visible(False)
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, ORANGE)]
    ax.legend(h, ['linear solves', 'assembly etc. (CPU ranks)'], ncol=2 if figsize[0] > 3 else 1,
              loc='lower left', bbox_to_anchor=(0, 1.0), handlelength=1.0, columnspacing=1.2)
    save(fig, name)


def breakdown():
    """One continental SSA stress balance (20 Picard its), unsmoothed start (H16_mps_it300)."""
    _breakdown([('MUMPS\n48 CPU cores (1 node)', 125, 12, '137 s  ·  4.2 SU'),
                ('hypre BoomerAMG\n1 A100 GPU + 16 cores', 66, 43, '109 s  ·  2.8 SU')],
               185, 25, 'breakdown.pdf')


def breakdown_transient():
    """Slide 11: one year of the SSA historical setup, 10 steps (TRM vs TRH_hostmpi). Linear solves
    include mass transport (MUMPS in both runs)."""
    # Compact: sits in a half-width column beside the slide's bullets.
    _breakdown([('MUMPS\n48 cores', 747, 84, '831 s\n23.0 SU'),
                ('BoomerAMG\nA100 + 16 cores', 276, 290, '565 s\n12.1 SU')],
               1050, 250, 'breakdown_transient.pdf', figsize=(2.35, 1.55))


def mps():
    """Slide 10: ranks sharing one A100 GPU, hypre BoomerAMG on the converged AIS3 matrix."""
    ranks = [1, 2, 4, 8, 16]
    default = [3.159, 13.844, 10.751, 15.100, 24.955]     # logs_181547397
    with_mps = [3.290, 5.359, 3.799, 4.048, 4.162]        # MPS=1 run
    fig, ax = plt.subplots(figsize=(2.55, 1.85))
    for y, c, lab in ((default, ORANGE, 'default'), (with_mps, BLUE, 'CUDA MPS')):
        ax.plot(ranks, y, color=c, linewidth=2, solid_capstyle='round', solid_joinstyle='round')
        ax.plot(ranks, y, 'o', color=c, markersize=5, markeredgecolor='white', markeredgewidth=1.2)
        ax.text(19, y[-1], f'{y[-1]:.1f} s', va='center', color=INK)
    ax.set_xscale('log', base=2)
    ax.set_xticks(ranks, [str(r) for r in ranks])
    ax.minorticks_off()
    ax.set_xlim(0.85, 34)
    ax.set_ylim(0, 27)
    ax.set_xlabel('MPI ranks sharing one A100 GPU')
    ax.set_ylabel('setup + solve (s)')
    grid(ax)
    h = [plt.Line2D([], [], color=c, linewidth=2, marker='o', markersize=4) for c in (ORANGE, BLUE)]
    ax.legend(h, ['default (time-sliced)', 'CUDA MPS (concurrent)'], loc='upper left', fontsize=8)
    save(fig, 'mps_ranks.pdf')


def iterations():
    """Slide 9: hypre BoomerAMG iterations per Picard iteration, 16-rank run (H16_mps_it300)."""
    its = {1: 200, 2: 300, 3: 300}
    its.update(zip(range(4, 21), [162, 60, 55, 51, 46, 46, 41, 35, 33, 30, 28, 23, 23, 21, 16, 14, 13]))
    failed = {1, 2}                                        # these handed over to MUMPS
    fig, ax = plt.subplots(figsize=(2.75, 1.75))
    for k, n in its.items():
        ax.bar(k, n, width=0.62, color=MUTED if k in failed else BLUE, edgecolor='white', linewidth=0.8)
    ax.axhline(300, color=AXIS, linewidth=0.8)
    ax.text(20.6, 300, 'cap 300', va='center', ha='right', fontsize=8, color=MUTED,
            bbox=dict(facecolor='white', edgecolor='none', pad=1))
    h = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, MUTED)]
    ax.legend(h, ['BoomerAMG solve', 'failed $\\rightarrow$ MUMPS,\nthen back'],
              loc='upper right', bbox_to_anchor=(1.0, 0.86), fontsize=8, handlelength=1.0)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.set_xlim(0.3, 20.7)
    ax.set_ylim(0, 330)
    ax.set_xlabel('Picard iteration')
    ax.set_ylabel('GMRES iterations')
    grid(ax)
    save(fig, 'picard_iterations.pdf')


def ais3_map():
    """Slide 7: observed speed, with the vertices behind the Picard-1 matrix's huge-diagonal cluster."""
    os.environ['ISSM_DIR'] = ('/g/data/vk83/apps/spack/1.1/release/linux-x86_64/'
                              'issm-git.2026.05.18_2026.05.18-kgta35igm37z4qnqnul7rcmgx2inftqd')
    import pyissm
    from matplotlib.tri import Triangulation
    dump = '/g/data/au88/jh7060/ACCESS-AIS3/execution_gpu_test/matrix_dump_ssa_picard1'
    with open(f'{dump}/A.bin', 'rb') as f:
        _, M, _, nz = np.fromfile(f, '>i4', 4)
        rl = np.fromfile(f, '>i4', M).astype(np.int64)
        cols = np.fromfile(f, '>i4', nz)
        vals = np.fromfile(f, '>f8', nz)
    rows = np.repeat(np.arange(M), rl)
    d = np.abs(vals[cols == rows])
    del cols, vals, rows
    v = np.load(f'{dump}/dof_vertex.npy')
    hot = np.unique(v[np.flatnonzero(d / np.median(d) > 1e5) // 2])
    md = pyissm.model.io.load_model('/g/data/au88/jh7060/ACCESS-AIS3/models/AIS3_SSA_inverted.nc')
    x = np.asarray(md.mesh.x).ravel() / 1e3
    y = np.asarray(md.mesh.y).ravel() / 1e3
    el = np.asarray(md.mesh.elements).astype(int) - 1
    vel = np.asarray(md.inversion.vel_obs, float).ravel()
    il = np.asarray(md.mask.ice_levelset).ravel()
    ol = np.asarray(md.mask.ocean_levelset).ravel()
    ice_el = (il[el] < 0).all(axis=1)
    tri = Triangulation(x, y, el[ice_el])
    cmap = LinearSegmentedColormap.from_list('seq_blue', SEQ)
    fig, ax = plt.subplots(figsize=(2.75, 2.55))
    pc = ax.tripcolor(tri, np.clip(vel, 1, 3000), cmap=cmap, norm=LogNorm(1, 3000),
                      shading='flat', rasterized=True)
    ax.tricontour(Triangulation(x, y, el), ol, levels=[0], colors=[INK2], linewidths=0.3)
    ax.plot(x[hot], y[hot], 'o', color=ORANGE, markersize=1.6, markeredgewidth=0, rasterized=True)
    ax.set_aspect('equal')
    ax.axis('off')
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.01)
    cb.set_label('observed speed (m/yr)', color=INK2, fontsize=8)
    cb.ax.tick_params(labelsize=8, color=AXIS)
    cb.outline.set_edgecolor(AXIS)
    ax.legend([plt.Line2D([], [], marker='o', linestyle='', color=ORANGE, markersize=4)],
              [f'{hot.size:,} vertices with\n$\\dot\\varepsilon\\approx0$ in Picard 1'],
              loc='lower left', fontsize=8, handletextpad=0.3, borderaxespad=0)
    print(f'hot vertices: {hot.size}')
    save(fig, 'ais3_zero_strain_map.pdf')


if __name__ == '__main__':
    charts = {'breakdown': breakdown, 'transient': breakdown_transient, 'mps': mps,
              'iterations': iterations}
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in charts:
        charts[what]()
    if what in ('charts', 'all'):
        for f in charts.values():
            f()
    if what in ('map', 'all'):
        ais3_map()
