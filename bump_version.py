#!/usr/bin/env python3

import semantic_version


def main():
    version = open('VERSION').read().strip()
    oldsemver = semantic_version.Version(version)
    newsemver = oldsemver.next_patch()
    print("Updating QuaTradis version from", oldsemver, "to", newsemver)
    open('VERSION', 'w').write(str(newsemver))


if __name__ == '__main__':
    main()