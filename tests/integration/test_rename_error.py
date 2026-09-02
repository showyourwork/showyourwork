import os

from helpers import TemporaryShowyourworkRepository

from showyourwork.config import edit_yaml

# A script that generates data and text file
# data file: read by plt script
# text file: read by variable command
data_gen_script = r"""
import paths
import numpy as np
import matplotlib.pyplot as plt

# data gen
N_pts=100
x=np.linspace(0, np.pi, N_pts)
y=np.sin(x)
data = np.column_stack((x, y))

# data save
np.savetxt(
    paths.data / "xy.dat",
    data,
    fmt="%.6e", # optional format
    header="x y", # optional header
    comments=""
)

# Write num points in a file
# will be called bz variable command
with open(paths.output / "points.txt", "w") as f:
    f.write(f"{N_pts}")
"""

# A script that plots generated data
fig_gen_script = r"""
import paths
import numpy as np
import matplotlib.pyplot as plt

# read generated data
x, y = np.loadtxt(
    paths.data / "xy.dat",
    unpack=True,
    skiprows=1,  # because of header="x y"
)

# plot and save
plt.plot(x,y)
plt.savefig(paths.figures / 'test_fig.png')
"""


class TestRenameError(TemporaryShowyourworkRepository):
    """Create the initial version of workflow:
    1. add python script for data, text gen
    2. add python script for fig gen
    3. edit tex to includegraphcs and variable
    5. edit Snakefile to add rule
    6. edit syw.yml to add dependency
    """

    def customize(self):
        """Create and edit all the necessary files for the workflow."""
        # Create the data_gen script
        with open(
            self.cwd / "src" / "scripts" / "data_gen.py", "w", encoding="utf-8"
        ) as f:
            print(data_gen_script, file=f)

        # Create the fig_gen script
        with open(
            self.cwd / "src" / "scripts" / "fig_gen.py", "w", encoding="utf-8"
        ) as f:
            print(fig_gen_script, file=f)

        # Edit tex to add variable command
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, encoding="utf-8") as f:
            ms_orig = f.read()
        with open(ms, "w", encoding="utf-8") as f:
            ms_new = ms_orig.replace(
                r"\end{document}",
                r"\begin{figure}"
                "\n"
                r"  \script{fig_gen.py}"
                "\n"
                r"  \includegraphics{figures/test_fig.png}"
                "\n"
                r"  \caption{Test figure}"
                "\n"
                r"  \label{fig:test_fig}"
                "\n"
                r"\end{figure}"
                "\n"
                r"The number of points used are "
                r"\variable{output/points.txt}."
                "\n"
                r"\end{document}",
            )
            print(ms_new, file=f)

        # Add a Snakemake rule
        sf = self.cwd / "Snakefile"
        with open(sf, encoding="utf-8") as f:
            contents = f.read()
        with open(sf, "w", encoding="utf-8") as f:
            print(contents, file=f)
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        "rule data_save:",
                        "    output:",
                        "        dat='src/data/xy.dat',",
                        "        out1='src/tex/output/points.txt'",
                        "    cache:",
                        "        True",
                        "    script:",
                        "        'src/scripts/data_gen.py'",
                    ]
                ),
                file=f,
            )

        # Edit showyourwork.yml to add dependency
        data_script_name = "src/scripts/data_gen.py"
        fig_script_name = "src/scripts/fig_gen.py"
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            if config.get("dependencies") is None:
                config["dependencies"] = {}
            config["dependencies"][fig_script_name] = [
                data_script_name,
                "src/data/xy.dat",
            ]

    def check_build(self):
        """
        1. Rename the txt file
            a. Snakefile
            b. script_file
            c. tex file
            d. delete old file
        2. rebuild
        """

        # Rename the txt file (points -> points_1)
        # # Edit Snakefile
        sf = self.cwd / "Snakefile"
        with open(sf, encoding="utf-8") as f:
            sf_orig = f.read()
        with open(sf, "w", encoding="utf-8") as f:
            sf_new = sf_orig.replace(
                r"points.txt",
                r"points_1.txt",
            )
            print(sf_new, file=f)

        # # Edit data gen script
        dat_gen_scrpt = self.cwd / "src" / "scripts" / "data_gen.py"
        with open(dat_gen_scrpt, encoding="utf-8") as f:
            dat_gen_scrpt_orig = f.read()
        with open(dat_gen_scrpt, "w", encoding="utf-8") as f:
            dat_gen_scrpt_new = dat_gen_scrpt_orig.replace(
                r"points.txt",
                r"points_1.txt",
            )
            print(dat_gen_scrpt_new, file=f)

        # # Edit tex
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, encoding="utf-8") as f:
            ms_orig = f.read()
        with open(ms, "w", encoding="utf-8") as f:
            ms_new = ms_orig.replace(
                r"points.txt",
                r"points_1.txt",
            )
            print(ms_new, file=f)

        # # Delete old file
        points_file = self.cwd / "src" / "tex" / "output" / "points.txt"
        if points_file.exists():
            os.remove(points_file)

        # Rebuild
        self.build_local()
