"""
Defines a custom exception for worfklow errors that
tries to print informative error messages at the very
end of the build log.

"""
import os
from pathlib import Path
import textwrap as tw

TEMPLATE = """
\033[0m{hline}

\033[1;{color}m{title}\033[0m

{hline}

\033[1;{color}mError:\033[0m {brief}
\033[1;{color}mRule:\033[0m {rule_name}
\033[1;{color}mScript:\033[0m {script}

\033[1;{color}mContext:\033[0m {context}

\033[1;{color}mDetailed message:\033[0m 

{message}
{hline}
"""


class ShowyourworkException(Exception):
    """
    Custom exception that stores error messages
    in a temporary file to be displayed at the
    end of the build process.

    """

    def __init__(
        self,
        message,
        exception_file=None,
        script=None,
        rule_name=None,
        context=None,
        delayed=True,
        brief="An error occurred while executing your workflow.",
        *args,
        **kwargs,
    ):
        # Record the inputs
        self.message = message
        self.script = script
        self.rule_name = rule_name
        self.context = context
        self.brief = brief

        # We only need to provide `exception_file` if we're not
        # inside the main Snakemake workflow run
        if exception_file is None:
            try:
                exception_file = files.exception
            except:
                exception_file = Path(".showyourwork/exception.log")

        # Format the info
        if script is None:
            script = "N/A"
        else:
            script = f"`showyourwork/workflow/scripts/{script}`"
        if rule_name is None:
            rule_name = "N/A"
        else:
            rule_name = f"{rule_name} in `showyourwork/workflow/rules/{rule_name}.smk`"
        if context is None:
            context = "N/A"
        CI = os.getenv("CI", "false") == "true"
        try:
            if not CI:
                width = os.get_terminal_size().columns
            else:
                width = 80
        except:
            width = 80
        if CI:
            color = "37"  # white
        else:
            color = "30"  # black
        hline = "*" * width
        title = "SHOWYOURWORK ERROR"
        pad = " " * max(0, (width - len(title)) // 2 - 2)
        title = f"{pad}{title}{pad}"
        full_message = TEMPLATE.format(
            hline=hline,
            title=title,
            script=script,
            rule_name=rule_name,
            brief="\n".join(
                tw.wrap(str(brief), width=width, initial_indent="       ")
            ).strip(),
            context="\n".join(
                tw.wrap(str(context), width=width, initial_indent="         ")
            ).strip(),
            message="\n".join(tw.wrap(str(message), width=width)),
            color=color,
        )

        if delayed:

            # Store the message in a temp file (read in the `onerror:` section
            # of the `Snakefile`, to be printed at the end of the log).
            with open(exception_file, "w") as f:
                print(full_message, file=f)
            super().__init__("\n\n" + full_message, *args, **kwargs)

        else:

            # Print it to the terminal right away
            print(full_message)
            super().__init__("\n\n" + message, *args, **kwargs)

    @staticmethod
    def print(exception_file=None):

        # We only need to provide `exception_file` if we're not
        # inside the main Snakemake workflow run
        if exception_file is None:
            try:
                exception_file = files.exception
            except:
                exception_file = Path(".showyourwork/exception.log")

        # Print any existing exceptions
        if exception_file.exists():
            with open(exception_file, "r") as f:
                message = f.read()
            print(message)
            exception_file.unlink()


class ShowyourworkWarning:
    """
    Custom warning that stores warning messages
    in a temporary file to be displayed at the
    end of the build process.

    """

    def __init__(
        self,
        message,
        exception_file=None,
        script=None,
        rule_name=None,
        context=None,
        brief="An error occurred while executing your workflow.",
    ):
        # Record the inputs
        self.message = message
        self.script = script
        self.rule_name = rule_name
        self.context = context
        self.brief = brief

        # We only need to provide `exception_file` if we're not
        # inside the main Snakemake workflow run
        if exception_file is None:
            try:
                exception_file = files.warning
            except:
                exception_file = Path(".showyourwork/warning.log")

        # Format the info
        if script is None:
            script = "N/A"
        else:
            script = f"`showyourwork/workflow/scripts/{script}`"
        if rule_name is None:
            rule_name = "N/A"
        else:
            rule_name = f"{rule_name} in `showyourwork/workflow/rules/{rule_name}.smk`"
        if context is None:
            context = "N/A"
        CI = os.getenv("CI", "false") == "true"
        try:
            if not CI:
                width = os.get_terminal_size().columns
            else:
                width = 80
        except:
            width = 80
        if CI:
            color = "37"  # white
        else:
            color = "30"  # black
        hline = "*" * width
        title = "SHOWYOURWORK WARNING"
        pad = " " * max(0, (width - len(title)) // 2 - 2)
        title = f"{pad}{title}{pad}"
        full_message = TEMPLATE.format(
            hline=hline,
            title=title,
            script=script,
            rule_name=rule_name,
            brief="\n".join(
                tw.wrap(str(brief), width=width, initial_indent="       ")
            ).strip(),
            context="\n".join(
                tw.wrap(str(context), width=width, initial_indent="         ")
            ).strip(),
            message="\n".join(tw.wrap(str(message), width=width)),
            color=color,
        )

        # Store the message in a temp file (read in the `onsuccess:` section
        # of the `Snakefile`, to be printed at the end of the log).
        with open(exception_file, "a+") as f:
            print(full_message, file=f)

    @staticmethod
    def print(exception_file=None):

        # We only need to provide `exception_file` if we're not
        # inside the main Snakemake workflow run
        if exception_file is None:
            try:
                exception_file = files.warning
            except:
                exception_file = Path(".showyourwork/warning.log")

        # Print any existing exceptions
        if exception_file.exists():
            with open(exception_file, "r") as f:
                message = f.read()
            print(message)
            exception_file.unlink()