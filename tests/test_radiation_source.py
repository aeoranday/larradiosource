from larradiosource.radiation import Source

from tomllib import load


def main():
    with open("../../config/radiation.toml", 'rb') as f:
        config = load(f)

    source = Source.model_validate(config['source'])
    decay_branch = source.get_random_decay_branch()
    print(decay_branch)
    return


if __name__ == "__main__":
    main()
