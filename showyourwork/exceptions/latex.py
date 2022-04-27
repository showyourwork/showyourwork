from .base import ShowyourworkException


class LaTeXException(ShowyourworkException):
    pass


class UnableToInferClassName(LaTeXException):
    def __init__(self, ms_name="ms"):
        super().__init__(
            f"Unable to determine document class in `{ms_name}.tex`."
        )


class TectonicError(LaTeXException):
    def __init__(self, logfile=None):
        if logfile:
            with open(logfile, "r") as f:
                tectonic_log = f.readlines()

            # Ensure the user imported showyourwork
            for line in tectonic_log:
                if "Package: showyourwork" in line:
                    showyourwork_imported = True
                    break
            else:
                showyourwork_imported = False

            if showyourwork_imported:

                # Scan the log for an error message
                for i, line in enumerate(tectonic_log[::-1]):
                    if line.startswith("!"):
                        message = "".join(tectonic_log[-i - 1 :])
                        break
                else:
                    message = (
                        "An error occurred while compiling the manuscript."
                    )
                message += f"\nFor more information, check out the log file:\n{logfile}."
            else:

                # Admonish the user (:
                message = r"Failed to compile manuscript. Did you forget to `\usepackage{showyourwork}`?"

        else:

            # No log to scan, so we're stuck with an uninformative message...
            message = "An error occurred while compiling the manuscript."

        super().__init__(message)


class FigureFormatError(LaTeXException):
    pass


class MissingXMLFile(LaTeXException):
    pass


class GraphicsPathError(LaTeXException):
    pass