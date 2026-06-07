import io
import sys

from config import RNACliqueConfig

from ruamel.yaml import YAML

def deep_dict_to_CommentedMap(yaml, d):
    # This might be slow, but it's concise and reliable.
    strio = io.StringIO()
    yaml.dump(d, strio)
    strio.seek(0)
    return yaml.load(strio)

def main():
    yaml = YAML()
    config = deep_dict_to_CommentedMap(yaml, RNACliqueConfig().marshal())
    for name, field in RNACliqueConfig.__dataclass_fields__.items():
        config.yaml_set_comment_before_after_key(
            name,
            field.metadata["description"]
        )
    yaml.dump(config, sys.stdout)
    

if __name__ == "__main__":
    main()
