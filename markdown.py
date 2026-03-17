import sys
import re
import textwrap
import contextlib

par_special_characters = re.compile(r"([#\[\]\\*_~^])")

class MarkdownDocument:
    def __init__(self, depth=1, file=sys.stdout, wrap=80):
        self._depth = depth
        self._file = file
        self._need_title = False
        self._wrapper = textwrap.TextWrapper(width=wrap)

    def enter_section(self, title):
        self._file.write("{} {}\n\n".format("#"*self._depth, title))
        self._need_title = False
        self._depth += 1

    def exit_section(self):
        self._depth -= 1
        self._need_title = True

    @classmethod
    def md_code(cls, s):
        assert isinstance(s, str)
        return "`{}`".format(s.replace("`", r"\`"))

    @contextlib.contextmanager
    def section(self, title):
        self.enter_section(title)
        try:
            yield self
        finally:
            self.exit_section()

    def write(self, contents, wrap=True):
        if self._need_title:
            self._file.write("{}\n\n".format("{}"*self._wrapper.width))
        if wrap:
            contents = self._wrapper.fill(contents)
        self._file.write(contents)
        
    def paragraph(self, contents, **kwargs):
        self.write(contents, **kwargs)
        self._file.write("\n\n")

    def code_block(self, contents, language=""):
        self.paragraph(f"```{language}\n{contents}\n```", wrap=False)

    @classmethod
    def escape(cls, text):
        return par_special_characters.sub(r"\\\1", text)
