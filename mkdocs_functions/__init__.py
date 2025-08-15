import os
from mkdocs_macros import fix_url

def define_env(env):
    @env.macro
    def doc_link(link, name=None, fix=True):
        if not name:
            name = link
        if fix:
            link = fix_url(link)
        return f"[{name}]({link})"

    @env.filter
    def comment_surround(x):
        return f"""-->
{x}
 <!--"""

    @env.macro
    def empty(x):
        return ""

    @env.filter
    def code_fence(x, language="text"):
        return f"""```{language}
{x}
```
"""
    @env.macro
    def git_branch():
        return os.environ.get("GITHUB_REF_NAME", "main")
    
    @env.macro
    def clone_command(branch):
        return (f"git clone --recurse-submodules -b {branch}"
                " https://github.com/actapia/rna_clique")
