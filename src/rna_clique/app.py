import sys
import functools
import os
import platform

from contextlib import contextmanager
from pathlib import Path

from . import config as config_module

eprint = functools.partial(print, file=sys.stderr)

class UnrecognizedFileExtensionError(Exception):
    def __init__(self, message, extension):
        super().__init__(message)
        self.extension = extension

def get_format_from_extension(extension_to_format):
    def get(export_out):
        extension = export_out.suffix[1:]
        try:
            return extension_to_format[extension]
        except KeyError:
            raise UnrecognizedFileExtensionError(
                f"Unrecognized file extension {extension!r}. "
                "Could not determine format.",
                extension
            )
    return get

def validate_input_dirs(config: config_module.RNACliqueConfig):
    try:
        config.validate_input_dirs()
    except config_module.InputValidationError:
        #eprint(f"{type(e).__name__}: {e}\n")
        eprint(
            ("There was a problem verifying the structure of the provided "
             "input.\n\nPlease check that each of your transcriptomes is in a "
             "separate directory named after the sample, and that all "
             "transcriptomes are named {}.\n").format(config.transcripts_name)
        )
        raise

# def validate_dir(path, name):
#     try:
#         config_module.RNACliqueConfig.validate_dir("
#x        )

transcript_id_parse_error_message = (
    "RNA-clique was unable to parse some of the FASTA headers in your data "
    "using the regular expression {}. Please check that all transcriptome "
    "transcript IDs can be parsed using the provided regular expression."
)

def print_transcript_id_parse_error_message(regex):
    eprint(
        f"{transcript_id_parse_error_message}\n".format(
            regex.pattern
        )
    )

def print_too_many_files_error_message(e):
    recommended_files = e.tried * 2            
    eprint("RNA-clique failed because it tried to open too many files.\n\n"
           f"Check that you can open at least {recommended_files}.\n")
    if platform.uname().system == "Darwin":
        eprint("If you are on macOS, you may need to raise the global "
               "per-process open file limit with: \n\n"
               "\tsudo sysctl -w kern.maxfilesperproc={0}\n\n"
               "If {0} also exceeds kern.maxfiles, you will also need to "
               "increase that value.\n".format(recommended_files))
    

generic_error_message = "An error occurred. {} did not finish successfully."

def print_generic_error_message(program):
    eprint(f"{generic_error_message}\n".format(program))

def except_hook(type_, value, traceback):
    print(f"{type_.__name__}: {value}", file=sys.stderr)

original_except_hook = sys.excepthook

def message_except_hook(message, original):
    def hook(*args, **kwargs):
        print(message, file=sys.stderr)
        original(*args, **kwargs)
    return hook



@contextmanager
def set_except_hook(
        verbose=bool(os.environ.get("RNA_CLIQUE_VERBOSE")),
        message=f"{generic_error_message}\n".format(Path(sys.argv[0]).name)
):
    previous_except_hook = sys.excepthook
    if verbose:
        sys.excepthook = message_except_hook(message, original_except_hook)
    else:
        sys.excepthook = message_except_hook(message, except_hook)
    yield
    sys.excepthook = previous_except_hook
