from helpers import TemporaryShowyourworkRepository
from showyourwork.config import edit_yaml

# A script that generates data and text file
# data file: read by plt script
# text file: read by variable command
data_gen_script = r"""
import paths
import numpy as np
import matplotlib.pyplot as plt

# Compute the age of the universe
N_pts=100
x=np.linspace(0, np.pi, N_pts)
y=np.sin(x)

# Combine x and y into a 2-column array
data = np.column_stack((x, y))

# plot and save
np.savetxt(
    paths.data / "xy.dat",
    data,
    fmt="%.6e", # optional format
    header="x y", # optional header
    comments=""
)

# Write num points in disk
with open(paths.output / "points.txt", "w") as f:
    f.write(f"{N_pts}")
"""

fig_gen_script = r"""
import paths
import numpy as np
import matplotlib.pyplot as plt

# Compute the age of the universe
x, y = np.loadtxt(
    paths.data / "xy.dat",
    unpack=True,
    skiprows=1,  # because of header="x y"
)

# plot and save
plt.plot(x,y)
plt.savefig(paths.figures / 'test_fig.png')
"""


class TestSnakeRule(TemporaryShowyourworkRepository):
    """Test the following:
    1. add python script for data, text gen 
    2. add python script for fig gen
    2. python script create cached data
    3. python script create output text
    4. python script plots cached data
    5. python script appends output text
    """

    # Keep this local only; remote covered elsewhere
    local_build_only = True

    def customize(self):
        """Create and edit all the necessary files for the workflow."""
        # Create the data_gen script
        with open(
            self.cwd / "src" / "scripts" / "data_gen.py",
              "w", encoding="utf-8"
              ) as f:
            print(data_gen_script, file=f)
        
        # Create the fig_gen script
        with open(
            self.cwd / "src" / "scripts" / "fig_gen.py",
              "w", encoding="utf-8"
              ) as f:
            print(fig_gen_script, file=f)

        # Import the variable into the tex file
        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"\end{document}",
                r"\begin{figure}" "\n"
                r"  \script{fig_gen.py}" "\n"
                r"  \includegraphics{figures/test_fig.png}" "\n"
                r"  \label{fig:test_fig}" "\n"
                r"\end{figure}" "\n"
                r"The number of points used are "
                r"\variable{output/points.txt}."
                "\n"
                r"\end{document}",
            )
            print(ms_new, file=f)

        # Add a Snakemake rule to run the script
        with open(self.cwd / "Snakefile", "r") as f:
            contents = f.read()
        with open(self.cwd / "Snakefile", "w") as f:
            print(contents, file=f)
            print("\n", file=f)
            print(
                "\n".join(
                    [
                        "rule data_save:",
                        "    output:",
                        "        dat='src/data/xy.dat',",
                        "        out1='src/tex/output/points.txt'",
                        "    cache: True",
                        "    script:",
                        "        'src/scripts/data_gen.py'",
                    ]
                ),
                file=f,
            )

        # showyourwork.yml
        data_script_name='src/scripts/data_gen.py'
        fig_script_name='src/scripts/fig_gen.py'
        with edit_yaml(self.cwd / "showyourwork.yml") as config:
            if config.get("dependencies") is None:
                config["dependencies"] = {}
            config["dependencies"][fig_script_name] = [
                data_script_name,
                'src/data/xy.dat']


    def check_build(self):
        """write something
        """
        # something
        # Import the variable into the tex file
        ms = self.cwd / "Snakefile"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"points.txt",
                r"points_1.txt",
            )
            print(ms_new, file=f)

        ms = self.cwd / "src" / "scripts" / "data_gen.py"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"points.txt",
                r"points_1.txt",
            )
            print(ms_new, file=f)

        ms = self.cwd / "src" / "tex" / "ms.tex"
        with open(ms, "r") as f:
            ms_orig = f.read()
        with open(ms, "w") as f:
            ms_new = ms_orig.replace(
                r"points.txt",
                r"points_1.txt",
            )
            print(ms_new, file=f)

        self.build_local()

        # # range output file name in snake
        # figure = self.cwd / "src" / "tex" / "figures" / "test_figure.pdf"
        # assert figure.exists()

        # # Check that the figure is present on the remote
        # for _n in range(self.auth_retries):
        #     try:
        #         overleaf.pull_files(
        #             [figure],
        #             self.overleaf_id,
        #             path=self.cwd,
        #             error_if_missing=True,
        #         )
        #     except exceptions.OverleafRateLimitExceeded:
        #         get_logger().warn(
        #             "Overleaf authentication failed. "
        #             f"Re-trying in {self.auth_sleep} seconds..."
        #         )
        #         time.sleep(self.auth_sleep)
        #     else:
        #         break

        # # Check that an exception is raised if we try to overwrite a file
        # # with uncommitted changes
        # ms = self.cwd / "src" / "tex" / "ms.tex"
        # with open(ms) as f:
        #     ms_orig = f.read()
        # with open(ms, "w") as f:
        #     f.write(r"% dummy comment\n" + ms_orig)
        # for _n in range(self.auth_retries):
        #     try:
        #         overleaf.pull_files(
        #             [ms],
        #             self.overleaf_id,
        #             path=self.cwd,
        #             error_if_local_changes=True,
        #         )
        #     except exceptions.OverleafRateLimitExceeded:
        #         get_logger().warn(
        #             f"Overleaf authentication failed. "
        #             f"Re-trying in {self.auth_sleep} seconds..."
        #         )
        #         time.sleep(self.auth_sleep)
        #     except exceptions.OverleafError:
        #         break
        #     else:
        #         raise Exception("Failed to raise exception!")

        # # Commit the changes and check that the exception is still raised
        # get_stdout(
        #     f'git add -f {ms} && git commit -m "changing ms.tex locally"',
        #     cwd=self.cwd,
        #     shell=True,
        # )
        # for _n in range(self.auth_retries):
        #     try:
        #         overleaf.pull_files(
        #             [ms],
        #             self.overleaf_id,
        #             path=self.cwd,
        #             error_if_local_changes=True,
        #         )
        #     except exceptions.OverleafRateLimitExceeded:
        #         get_logger().warn(
        #             "Overleaf authentication failed. "
        #             f"Re-trying in {self.auth_sleep} seconds..."
        #         )
        #         time.sleep(self.auth_sleep)
        #     except exceptions.OverleafError:
        #         break
        #     else:
        #         raise Exception("Failed to raise exception!")

        # # Amend the commit message with the magical `[showyourwork]` label
        # # and check that the merge works
        # get_stdout(
        #     'git commit --amend -m "[showyourwork] changing ms.tex locally"',
        #     cwd=self.cwd,
        #     shell=True,
        # )
        # for _n in range(self.auth_retries):
        #     try:
        #         overleaf.pull_files(
        #             [ms],
        #             self.overleaf_id,
        #             path=self.cwd,
        #             error_if_local_changes=True,
        #         )
        #     except exceptions.OverleafRateLimitExceeded:
        #         get_logger().warn(
        #             "Overleaf authentication failed. "
        #             f"Re-trying in {self.auth_sleep} seconds..."
        #         )
        #         time.sleep(self.auth_sleep)
        #     else:
        #         break

        # # ---- Regression test for issue #603 --------------------------------
        # # After the first successful build, simulate a *new* Overleaf-side
        # # edit and rebuild WITHOUT cleaning.  Before #603 was fixed the
        # # Overleaf pull lived inside the ``onstart`` handler, which Snakemake
        # # skips when it considers everything up-to-date, so the changes were
        # # silently ignored.
        # ms_after_first_build = ms.read_text()
        # marker = "% issue-603-test-marker"
        # ms.write_text(ms_after_first_build + "\n" + marker + "\n")

        # for _n in range(self.auth_retries):
        #     try:
        #         overleaf.push_files(
        #             [ms],
        #             self.overleaf_id,
        #             path=self.cwd,
        #         )
        #     except exceptions.OverleafRateLimitExceeded:
        #         get_logger().warn(
        #             "Overleaf authentication failed. "
        #             f"Re-trying in {self.auth_sleep} seconds..."
        #         )
        #         time.sleep(self.auth_sleep)
        #     else:
        #         break

        # # Revert locally so the rebuild must pull from Overleaf
        # ms.write_text(ms_after_first_build)
        # get_stdout("git checkout -- src/tex/ms.tex", shell=True, cwd=self.cwd)

        # # Rebuild without cleaning
        # self.build_local()

        # # The marker must now be present in the local manuscript
        # assert (
        #     marker in ms.read_text()
        # ), "Second build did not pull the new Overleaf changes (issue #603)"

# class TestSnakeRule(TemporaryShowyourworkRepository):
#     local_build_only = False