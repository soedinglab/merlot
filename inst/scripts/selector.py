import argparse

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('script', help='Which script to use.')
    args = parser.parse_args()


if __name__ == "__main__":
    main()
